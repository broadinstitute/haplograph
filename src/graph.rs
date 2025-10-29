
use log::{info, warn};
use anyhow::{Result as AnyhowResult, Context};
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use std::path::{PathBuf};
use std::collections::{HashMap, BTreeMap, HashSet};
use std::io::Write;
use std::fs::File;
use crate::intervals;
use crate::util;
use crate::methyl;
use std::error::Error;
use indicatif::{ProgressBar, ProgressStyle};
use bio::io::fastq;


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

pub fn get_node_edge_info(windows: &Vec<(String, usize, usize)>, final_hap_list: &Vec<HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)>>, min_reads: usize, frequency_min: f64, haplotype_number: usize) -> (HashMap<String, NodeInfo>, HashMap<(String,String), EdgeInfo>) {
     let mut node_info = HashMap::new();
    let mut edge_info = HashMap::new();

    for (index, window) in windows.iter().enumerate() {
        let final_hap_list_index = final_hap_list[index].clone();

        for (i, (final_haplotype_seq,(cigar, read_dict, allele_frequency))) in final_hap_list_index.iter().enumerate(){
            let read_vector = read_dict.keys().cloned().collect::<Vec<_>>();
            let read_vector_clone = read_vector.iter().map(|x| x.split("|").collect::<Vec<_>>()[0].to_string()).collect::<HashSet<_>>();
            let node_id = format!("H.{}:{}-{}.{}", window.0, window.1, window.2, i);
            // let cigar = cigar_dict_list[index].get(final_haplotype_seq).unwrap();
            let read_vector_len = read_vector.len();
            node_info.insert(node_id.clone(), 
             NodeInfo {
                nodename: node_id.clone(),
                pos: window.1,
                seq: final_haplotype_seq.clone(),
                cigar: cigar.clone(),
                support_reads: read_vector_len,
                allele_frequency: allele_frequency.clone(),
                methyl_info: read_dict.clone(),
                haplotype_index: None,
            });
            
            // add edge information
            if index < windows.len() - 1 {
                let next_window = windows[index + 1].clone();
                // every node only have one choice for the next window
                for (j, (next_final_haplotype_seq, (next_cigar,next_methyl_dict, next_allele_frequency))) in final_hap_list[index + 1].iter().enumerate(){
                    let next_node_id = format!("H.{}:{}-{}.{}", next_window.0, next_window.1, next_window.2, j);
                    // let next_cigar = cigar_dict_list[index + 1].get(next_final_haplotype_seq).unwrap();
                    let next_read_vector = next_methyl_dict.keys().cloned().collect::<Vec<_>>();
                    let next_read_vector_len = next_read_vector.len();
                    let next_read_vector_clone = next_read_vector.iter().map(|x| x.split("|").collect::<Vec<_>>()[0].to_string()).collect::<HashSet<_>>();
                    let overlapping_reads = read_vector_clone.intersection(&next_read_vector_clone).cloned().collect::<Vec<_>>();
                    let overlap_ratio = overlapping_reads.len() as f64 / (read_vector_len as f64).max( next_read_vector_len as f64);
                    // println!("readset1: {}, readset2: {}, overlap_ratio: {}", read_vector.len(), next_read_vector.len(), overlap_ratio);
                    if overlapping_reads.len() >= min_reads - 1 {
                        edge_info.insert((node_id.clone(), next_node_id.clone()), 
                        EdgeInfo {
                            src: node_id.clone(),
                            dst: next_node_id.clone(),
                            overlap_ratio: overlap_ratio,
                            overlapping_reads: overlapping_reads.join(","),
                        });
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
        // println!("methylation_info, {:?}", methyl_info);
        let read_names = node_info.methyl_info.keys().cloned().collect::<Vec<_>>().join(",");
        let mod_score_dict = methyl::aggregate_methylation_reads(methyl_info, methyl_threshold);
        let mod_score_dict_string = mod_score_dict.iter().map(|(pos, score)| format!("{}:{}", pos, score)).collect::<Vec<_>>().join(",");
        let average_mod_score = mod_score_dict.values().sum::<f32>() / mod_score_dict.len() as f32;
        // println!("mod_score_dict: {:?}", mod_score_dict);

        let mut node_info_clone = BTreeMap::new();
        node_info_clone.insert("pos".to_string(), start.to_string());
        // node_info_clone.insert("seq".to_string(), haplotype_seq);
        node_info_clone.insert("cigar".to_string(), haplotype_cigar.to_string());
        node_info_clone.insert("support_reads".to_string(), read_num.to_string());
        node_info_clone.insert("allele_frequency".to_string(), allele_frequency.to_string());
        node_info_clone.insert("read_names".to_string(), read_names.to_string());
        node_info_clone.insert("mod_score_dict".to_string(), mod_score_dict_string);
        // node_info_clone.insert("read_names".to_string(), read_names);
        // anchor_info_clone.seq = String::new();
        let json_string =
            serde_json::to_string(&node_info_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = format!("S\t{}\t{}\tPG:J:{}\tRC:i:{}\tML:f:{:.2}", haplotype_id, haplotype_seq, json_string, read_num, average_mod_score);
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


pub fn start(bam: &mut IndexedReader, windows: &Vec<(String,usize, usize)>, reference_fa: &Vec<fastq::Record>, sampleid: &String, min_reads: usize, methyl_threshold: f32, frequency_min: f64, primary_only: bool, output_prefix: &String, haplotype_number: usize) -> AnyhowResult<()> {

    let mut final_hap_list = Vec::new(); 
    for (i, window) in windows.iter().enumerate() {
        let (chromosome, start, end) = window;
        let haplotype_info = intervals::start(bam, &reference_fa, &chromosome, *start, *end, &sampleid, min_reads as usize, frequency_min, primary_only, false).unwrap();
        final_hap_list.push(haplotype_info);
    }

    let (node_info, edge_info) = get_node_edge_info(&windows, &final_hap_list, 1 as usize, frequency_min, haplotype_number);
    
    let gfa_output = PathBuf::from(format!("{}.gfa", output_prefix));
    let _ = write_gfa_output(&node_info, &edge_info, &gfa_output, methyl_threshold);

    info!("Graph reconstruction completed");
    Ok(())
}