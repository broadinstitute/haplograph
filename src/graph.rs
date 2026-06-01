use crate::intervals;
use crate::methyl;
use crate::util;
use anyhow::{Context, Result as AnyhowResult};
use bio::io::fastq;
use log::info;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

pub struct NodeInfo {
    pub nodename: String,
    pub pos: usize,
    pub seq: String,
    pub cigar: String,
    pub support_reads: usize,
    pub allele_frequency: f64,
    pub methyl_info: HashMap<String, HashMap<usize, f32>>,
    pub haplotype_index: Option<Vec<usize>>,
}
pub struct EdgeInfo {
    pub src: String,
    pub dst: String,
    pub overlap_ratio: f64,
    pub overlapping_reads: String,
}

/// graph connectivity based on reference alignment order
pub fn get_node_edge_info(
    windows: &Vec<(String, usize, usize)>,
    final_hap_list: &Vec<(
        HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)>,
        HashMap<String, (u64, u64)>,
        HashMap<String, Vec<u8>>,
        HashMap<String, String>,
    )>,
    min_reads: usize
) -> (
    HashMap<String, NodeInfo>,
    HashMap<(String, String), EdgeInfo>,
) {
    let mut node_info = HashMap::new();
    let mut edge_info = HashMap::new();

    for (index, window) in windows.iter().enumerate() {
        let (
            final_hap_list_index,
            coordinate_list_index,
            _read_sequence_dict_index,
            _read_strand_dict_index,
        ) = final_hap_list[index].clone();

        for (i, (final_haplotype_seq, (cigar, read_dict, allele_frequency))) in
            final_hap_list_index.iter().enumerate()
        {
            let read_vector = read_dict.keys().cloned().collect::<Vec<_>>();
            let read_vector_clone = read_vector
                .iter()
                .map(|x| x.split("|").collect::<Vec<_>>()[0].to_string())
                .collect::<HashSet<_>>();
            let node_id = format!("H.{}:{}-{}.{}", window.0, window.1, window.2, i);
            let read_vector_len = read_vector.len();
            node_info.insert(
                node_id.clone(),
                NodeInfo {
                    nodename: node_id.clone(),
                    pos: window.1,
                    seq: final_haplotype_seq.clone(),
                    cigar: cigar.clone(),
                    support_reads: read_vector_len,
                    allele_frequency: *allele_frequency,
                    methyl_info: read_dict.clone(),
                    haplotype_index: None,
                },
            );

            if index < windows.len() - 1 {
                let next_window = windows[index + 1].clone();
                for (
                    j,
                    (
                        next_final_haplotype_seq,
                        (next_cigar, next_methyl_dict, next_allele_frequency),
                    ),
                ) in final_hap_list[index + 1].0.iter().enumerate()
                {
                    let next_node_id = format!(
                        "H.{}:{}-{}.{}",
                        next_window.0, next_window.1, next_window.2, j
                    );
                    let next_read_vector = next_methyl_dict.keys().cloned().collect::<Vec<_>>();
                    let next_read_vector_len = next_read_vector.len();
                    let next_read_vector_clone = next_read_vector
                        .iter()
                        .map(|x| x.split("|").collect::<Vec<_>>()[0].to_string())
                        .collect::<HashSet<_>>();
                    let overlapping_reads = read_vector_clone
                        .intersection(&next_read_vector_clone)
                        .cloned()
                        .collect::<Vec<_>>();
                    let overlap_ratio = overlapping_reads.len() as f64
                        / (read_vector_len as f64).max(next_read_vector_len as f64);
                    if overlapping_reads.len() >= min_reads - 1 {
                        edge_info.insert(
                            (node_id.clone(), next_node_id.clone()),
                            EdgeInfo {
                                src: node_id.clone(),
                                dst: next_node_id.clone(),
                                overlap_ratio,
                                overlapping_reads: overlapping_reads.join(","),
                            },
                        );
                    }
                }
            }
        }
    }
    (node_info, edge_info)
}

pub fn write_gfa_output(
    node_file: &HashMap<String, NodeInfo>,
    edge_info: &HashMap<(String, String), EdgeInfo>,
    output_filename: &PathBuf,
    methyl_threshold: f32,
) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;
    writeln!(file, "H\tVN:Z:1.0")?;

    let mut node_output = Vec::new();
    for (haplotype_id, node_info) in node_file.iter() {
        let hap_name = haplotype_id.split(".").nth(1).unwrap().to_string();
        let (chromosome, start, end) = util::split_locus(hap_name);
        let haplotype_seq = node_info.seq.clone();
        let haplotype_cigar = node_info.cigar.clone();
        let read_num = node_info.support_reads;
        let allele_frequency = node_info.allele_frequency;
        let methyl_info = node_info.methyl_info.clone();
        let read_names = node_info
            .methyl_info
            .keys()
            .cloned()
            .collect::<Vec<_>>()
            .join(",");
        let mod_score_dict = methyl::aggregate_methylation_reads(methyl_info, methyl_threshold);
        let mod_score_dict_string = mod_score_dict
            .iter()
            .map(|(pos, score)| format!("{}:{}", pos, score))
            .collect::<Vec<_>>()
            .join(",");
        let average_mod_score = if mod_score_dict.is_empty() {
            0.0
        } else {
            mod_score_dict.values().sum::<f32>() / mod_score_dict.len() as f32
        };
        let mut node_info_clone = BTreeMap::new();
        node_info_clone.insert("pos".to_string(), start.to_string());
        node_info_clone.insert("cigar".to_string(), haplotype_cigar.to_string());
        node_info_clone.insert("support_reads".to_string(), read_num.to_string());
        node_info_clone.insert(
            "allele_frequency".to_string(),
            format!("{:.4}", allele_frequency),
        );
        node_info_clone.insert("read_names".to_string(), read_names.to_string());
        node_info_clone.insert("mod_score_dict".to_string(), mod_score_dict_string);
        let json_string =
            serde_json::to_string(&node_info_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = format!(
            "S\t{}\t{}\tPG:J:{}\tRC:i:{}\tML:f:{:.2}",
            haplotype_id, haplotype_seq, json_string, read_num, average_mod_score
        );
        node_output.push(formatted_string);
    }

    let mut link_output = Vec::new();

    for ((src, dst), edge_info) in edge_info.iter() {
        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", src, dst));
    }
    node_output.sort();
    link_output.sort();

    for s in node_output {
        writeln!(file, "{}", s)?;
    }
    for s in link_output {
        writeln!(file, "{}", s)?;
    }

    Ok(())
}

pub fn start(
    bam_path: &String,
    windows: &Vec<(String, usize, usize)>,
    reference_fa: &Vec<fastq::Record>,
    sampleid: &String,
    min_reads: usize,
    methyl_threshold: f32,
    frequency_min: f64,
    primary_only: bool,
    pileup:bool,
    output_prefix: &String,
    minimal_gap_length: usize,
) -> AnyhowResult<()> {
    // Parallelize window processing - each thread gets its own BAM reader
    info!("Processing {} windows in parallel", windows.len());
    let final_hap_list: Vec<_> = windows
        .par_iter()
        .map(|window| {
            let (chromosome, start, end) = window;
            // Create a new BAM reader for this thread
            let mut bam = util::open_bam_file(bam_path);
            // Process this window
            intervals::start(
                &mut bam,
                reference_fa,
                chromosome,
                *start,
                *end,
                sampleid,
                min_reads,
                frequency_min,
                primary_only,
                pileup,
                false
            )
            .with_context(|| format!("Failed to process window {}:{}-{}", chromosome, start, end))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let (node_info, edge_info) = get_node_edge_info(windows, &final_hap_list, 1_usize);
    println!("node_info: {:?}", node_info.len());
    println!("edge_info: {:?}", edge_info.len());
    let gfa_output = PathBuf::from(format!("{}.gfa", output_prefix));
    let _ = write_gfa_output(&node_info, &edge_info, &gfa_output, methyl_threshold);

    info!("Graph reconstruction completed");
    Ok(())
}

#[cfg(test)]
mod graph_e2e_tests {
    use super::start as graph_start;
    use bio::io::fastq;
    use std::path::Path;

    #[test]
    fn graph_build_produces_gfa_nodes_on_tp53_window() {
        let bam_path = "HG00097_tp53_test.bam";
        if !Path::new(bam_path).exists() {
            return;
        }

        let reference = vec![fastq::Record::with_attrs(
            "chr17",
            None,
            vec![b'N'; 10_000_000].as_slice(),
            &[],
        )];
        let windows = vec![("chr17".to_string(), 7668752_usize, 7668828_usize)];
        let output_prefix = "test_tp53_graph_e2e";

        graph_start(
            &bam_path.to_string(),
            &windows,
            &reference,
            &"HG00097".to_string(),
            2,
            0.5,
            0.0,
            false,
            true,
            &output_prefix.to_string(),
            1,
        )
        .expect("graph build should succeed");

        let gfa_path = format!("{output_prefix}.gfa");
        let gfa = std::fs::read_to_string(&gfa_path).expect("GFA should be written");
        let s_line_count = gfa.lines().filter(|line| line.starts_with('S')).count();

        let _ = std::fs::remove_file(&gfa_path);

        assert!(
            s_line_count > 0,
            "expected GFA segment lines after pileup sequence fix, got {s_line_count}"
        );
    }
}
