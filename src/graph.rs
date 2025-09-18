
use log::{info, warn};
use anyhow::{Result as AnyhowResult, Context};
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use std::path::{PathBuf};
use std::collections::{HashMap, BTreeMap, HashSet};
use std::io::Write;
use std::fs::File;
use crate::intervals;
use crate::util;
use std::error::Error;
use indicatif::{ProgressBar, ProgressStyle};
use bio::io::fastq;


pub fn get_node_edge_info(locus: &String, windows: &Vec<(usize, usize)>, final_hap_list: &Vec<HashMap<String, (String,Vec<String>, f64)>>, min_reads: usize, frequency_min: f64) -> (HashMap<String, String>, HashMap<(String,String), String>) {
    let (chromosome, start, end) = util::split_locus(locus.clone());
     let mut node_info = HashMap::new();
    let mut edge_info = HashMap::new();
    for (index, window) in windows.iter().enumerate() {
        for (i, (final_haplotype_seq,(cigar, read_vector, allele_frequency))) in final_hap_list[index].iter().enumerate(){
            let haplotype_id = format!("H.{}:{}-{}.{}", chromosome, window.0, window.1, i);
            // let cigar = cigar_dict_list[index].get(final_haplotype_seq).unwrap();
            let read_vector_len = read_vector.len();
            node_info.insert(haplotype_id.clone(), format!("L\t{}\t{}\t{}\t{:.2}\t{}", final_haplotype_seq, cigar, read_vector_len, allele_frequency, read_vector.join(",")));

            // add edge information
            if index < windows.len() - 1 {
                let next_window = windows[index + 1];
                for (j, (next_final_haplotype_seq, (next_cigar,next_read_vector, next_allele_frequency))) in final_hap_list[index + 1].iter().enumerate(){
                    let next_haplotype_id = format!("H.{}:{}-{}.{}", chromosome, next_window.0, next_window.1, j);
                    // let next_cigar = cigar_dict_list[index + 1].get(next_final_haplotype_seq).unwrap();
                    let next_read_vector_len = next_read_vector.len();
                    let overlapping_reads = util::find_overlapping_reads(read_vector, next_read_vector);
                    let overlap_ratio = overlapping_reads.len() as f64 / (read_vector_len as f64).min( next_read_vector_len as f64);
                    // println!("readset1: {}, readset2: {}, overlap_ratio: {}", read_vector.len(), next_read_vector.len(), overlap_ratio);
                    if overlapping_reads.len() >= 1 {
                        edge_info.insert((haplotype_id.clone(), next_haplotype_id.clone()), format!("E\t{}\t{}\tOverlapRatio:{:.3}\tOverlappingReads:{}", 
                            read_vector_len, next_read_vector_len, overlap_ratio, overlapping_reads.join(",")));
                    }
                    
                }
            }
        }
    }
    (node_info, edge_info)
}

pub fn write_gfa_output(
    node_file: &HashMap<String, String>,
    edge_info: &HashMap<(String, String), String>,
    output_filename: &PathBuf,
) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;
    writeln!(file, "H\tVN:Z:1.0")?;

    let mut node_output = Vec::new();
    for (haplotype_id, node_info) in node_file.iter() {
        let hap_name = haplotype_id.split("|").next().unwrap().to_string();
        let pos_string = hap_name.split(".").next().unwrap().to_string();
        let (chromosome, start, end) = util::split_locus(pos_string);
        let parts: Vec<&str> = node_info.split("\t").collect();
        if parts.len() >= 5 {
            let haplotype_seq = parts[1];
            let haplotype_cigar = parts[2];
            let read_num = parts[3];
            let allele_frequency = parts[4];
            
            let mut node_info_clone = BTreeMap::new();
            node_info_clone.insert("pos".to_string(), start.to_string());
            // node_info_clone.insert("seq".to_string(), haplotype_seq);
            node_info_clone.insert("cigar".to_string(), haplotype_cigar.to_string());
            node_info_clone.insert("support_reads".to_string(), read_num.to_string());
            node_info_clone.insert("allele_frequency".to_string(), allele_frequency.to_string());
            // node_info_clone.insert("read_names".to_string(), read_names);
            // anchor_info_clone.seq = String::new();
            let json_string =
                serde_json::to_string(&node_info_clone).unwrap_or_else(|_| "{}".to_string());
            let formatted_string = format!("S\t{}\t{}\tPG:J:{}", haplotype_id, haplotype_seq, json_string);
            node_output.push(formatted_string);
        }

    }

    let mut link_output = Vec::new();

    for ((src, dst), edge_info_string) in edge_info.iter() {

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


pub fn start(bam: &mut IndexedReader, windows: &Vec<(usize, usize)>, locus: &String, reference_fa: &Vec<fastq::Record>, sampleid: &String, min_reads: usize, frequency_min: f64, primary_only: bool, output_prefix: &String) -> AnyhowResult<()> {

    let (chromosome, start, end) = util::split_locus(locus.clone());
    let mut final_hap_list = Vec::new(); 
    for (i, window) in windows.iter().enumerate() {
        let (start, end) = window;
        let haplotype_info = intervals::start(bam, &reference_fa, &chromosome, *start, *end, &sampleid, min_reads as usize, frequency_min, primary_only, false).unwrap();
        final_hap_list.push(haplotype_info);
    }

    let (node_info, edge_info) = get_node_edge_info(&locus, &windows, &final_hap_list, min_reads as usize, frequency_min);
    // let gfa_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplograph.gfa", cli.output, cli.sampleid, chromosome, start, end));
 
    let gfa_output = PathBuf::from(format!("{}.gfa", output_prefix));
    write_gfa_output(&node_info, &edge_info, &gfa_output);

    info!("Graph reconstruction completed");
    Ok(())
}