use log::{info, warn};
use anyhow::{Result as AnyhowResult, Context};
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use std::path::{PathBuf};
use std::collections::{HashMap, BTreeMap, HashSet};
use std::io::Write;
use std::fs::File;
use std::io::{BufRead, BufReader};
use flate2::read::GzDecoder;
use serde_json::Value;
use std::error::Error;
use crate::eval;
use crate::util;
pub struct NodeInfo {
    pub seq: String,
    pub cigar: String,
    pub support_reads: String,
    pub allele_frequency: String,
    pub read_names: String,
}

pub fn load_graph(filename: &PathBuf) -> AnyhowResult< (HashMap<String, NodeInfo>,  HashMap<String, Vec<String>>) > {
    let file = File::open(filename)?;
    let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut node_info = HashMap::new();
    let mut edge_info = HashMap::new();

    for line in reader.lines() {
        let line = line?.trim_end().to_string();
        if line.starts_with('S') {
            let itemlist: Vec<&str> = line.split('\t').collect();
            let name = itemlist[1];
            let seq = itemlist[2];
            let annotation = &itemlist[3][5..];

            let json_value: Value = serde_json::from_str(annotation).unwrap();

            let mut value = json_value;

            let value_info = NodeInfo {
                seq: seq.to_string(),
                cigar: value["cigar"].to_string(),
                support_reads: value["support_reads"].to_string(),
                allele_frequency: value["allele_frequency"].to_string(),
                read_names: value["read_names"].to_string(),
            };
            node_info.insert(name.to_string(), value_info);

        } else if line.starts_with('L') {
            let itemlist: Vec<&str> = line.split('\t').collect();
            let src = itemlist[1];
            let dst = itemlist[3];
            edge_info.entry(src.to_string()).or_insert(Vec::new()).push(dst.to_string());
        }
    }
    Ok((node_info, edge_info))
}

pub fn find_source_node(edge_info: &HashMap<String, Vec<String>>) -> Vec<String> {
    let mut source_nodes = HashSet::new();
    let mut target_nodes = HashSet::new();
    for (src, dst_list) in edge_info.iter() {
        source_nodes.insert(src.clone());
        for dst in dst_list.iter() {
            target_nodes.insert(dst.clone());
        }
    }
    let start_nodes = source_nodes.difference(&target_nodes).collect::<HashSet<&String>>();
    start_nodes.into_iter().map(|s| s.clone()).collect()
}

/// Enumerate all paths in a directed acyclic graph using DFS
pub fn enumerate_all_paths(
    _node_info: &HashMap<String, NodeInfo>, 
    edge_info: &HashMap<String, Vec<String>>
) -> std::result::Result<Vec<Vec<String>>, Box<dyn Error>> {
    let source_nodes = find_source_node(edge_info);
    let mut all_paths = Vec::new();
    
    for src in source_nodes {
        let mut path = Vec::new();
        path.push(src.clone());
        dfs_traverse(&src, edge_info, &mut path, &mut all_paths);
    }
    
    Ok(all_paths)
}

/// Recursive DFS to find all paths from a starting node
fn dfs_traverse(
    current_node: &String,
    edge_info: &HashMap<String, Vec<String>>,
    current_path: &mut Vec<String>,
    all_paths: &mut Vec<Vec<String>>
) {
    // If no outgoing edges, this is a complete path
    if !edge_info.contains_key(current_node) {
        all_paths.push(current_path.clone());
        return;
    }
    
    // Explore all outgoing edges
    if let Some(next_nodes) = edge_info.get(current_node) {
        for next_node in next_nodes {
            current_path.push(next_node.clone());
            dfs_traverse(next_node, edge_info, current_path, all_paths);
            current_path.pop(); // Backtrack
        }
    }
}

pub fn traverse_graph(node_info: &HashMap<String, NodeInfo>, edge_info: &HashMap<String, Vec<String>>) -> std::result::Result<HashMap<Vec<String>, (String, HashSet<String>)>, Box<dyn Error>> {
    let all_paths = enumerate_all_paths(node_info, edge_info)?;
    
    // Get the sequence of each path
    let mut all_sequences = HashMap::new();
    for path in all_paths {
        let mut sequence = String::new();
        // let mut total_supported_reads = 0;
        let mut read_names = HashSet::new();
        for node in path.iter() {
            let node_info_dict = node_info.get(node).unwrap();
            // let supported_reads = node_info_dict.support_reads.clone();
            // let supported_reads = supported_reads.trim_matches('"').parse::<usize>().unwrap();
            let supported_reads = node_info_dict.read_names.clone();
            for read in supported_reads.split(",").collect::<Vec<_>>(){
                read_names.insert(read.to_string());
            }
            let haplotype_seq = node_info_dict.seq.clone();
            sequence += &haplotype_seq;
            // total_supported_reads += supported_reads;
        }
        all_sequences.insert(path.clone(), (sequence, read_names.clone()));
    }

    Ok(all_sequences)
}

pub fn write_graph_path_fasta(all_sequences: &HashMap<Vec<String>, (String, HashSet<String>)>, output_filename: &PathBuf) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;
    let chars_per_line = 60;
    for (index,(path, (sequence, supported_reads))) in all_sequences.iter().enumerate() {
        let chromosome = path[0].split(".").collect::<Vec<_>>()[1].split(":").collect::<Vec<_>>()[0];
        let (start, end) = eval::find_alignment_intervals(path.iter().map(|x| x.as_str()).collect::<Vec<_>>())?;
        writeln!(file, ">{}:{}-{}.{}\tSupports:{}\t{}", chromosome, start, end, index, supported_reads.len(), path.join("|"))?;
        // write the sequence in fasta format
        let seq_len = sequence.len();
        let full_lines = seq_len / chars_per_line;
        for i in 0..full_lines {
            let start = i * chars_per_line;
            let end = start + chars_per_line;
            writeln!(file, "{}", &sequence[start..end])?;
        }
        // Write any remaining characters that didn't make up a full line
        if seq_len % chars_per_line != 0 {
            writeln!(file, "{}", &sequence[full_lines * chars_per_line..])?;
        }
    }
    Ok(())
}

pub fn start(graph_filename: &PathBuf, locus: &String, germline_only:bool, haplotype_number: usize, output_prefix: &PathBuf) -> AnyhowResult<()> {
    let (node_info, edge_info) = load_graph(graph_filename).unwrap();
    let all_sequences = traverse_graph(&node_info, &edge_info).unwrap();
    let (ref_chromosome, ref_start, ref_end) = util::split_locus(locus.to_string());
    let mut hap_number = haplotype_number;
    let all_sequences_filtered = if germline_only {
        // sort all_sequences by the spanning length and the supported_reads
        let mut full_sequences = Vec::new();
        let mut all_reads = HashSet::new();
        for (path, (sequence, supported_reads)) in all_sequences.iter() {
            all_reads.extend(supported_reads.iter().cloned());
            let (start, end) = eval::find_alignment_intervals(path.iter().map(|x| x.as_str()).collect::<Vec<_>>())?;
            if start == ref_start && end == ref_end {
                full_sequences.push((path.clone(), (sequence.clone(), supported_reads.clone())));
            }
        }
        full_sequences.sort_by(|a, b| b.1.1.len().cmp(&a.1.1.len()));
        if &full_sequences.len() < &haplotype_number {
            hap_number = full_sequences.len();
        }
        // find the first haplotypes
        let primary_haplotype_1: HashMap<Vec<String>, (String, HashSet<String>)> = full_sequences.clone().into_iter().take(1).map(|(k, v)| (k.clone(), v.clone())).collect();
        let read_set_hap1 = primary_haplotype_1.values().map(|(_, supported_reads)| supported_reads).collect::<Vec<_>>()[0];
        let mut rest_primary_haplotypes: Vec<(Vec<String>, (String, HashSet<String>))> = Vec::new();
        let read_set_hap2: HashSet<String> = (&all_reads - read_set_hap1).iter().cloned().collect();
        for (path, (sequence, supported_reads)) in full_sequences.into_iter().skip(1) {
            let overlapped_reads: HashSet<String> = supported_reads.intersection(&read_set_hap2).cloned().collect();
            rest_primary_haplotypes.push((path.clone(), (sequence.clone(), overlapped_reads)));
        }
        // sort the primary_haplotypes by the supported_reads
        rest_primary_haplotypes.sort_by(|a, b| b.1.1.len().cmp(&a.1.1.len()));
        let primary_haplotypes_2: HashMap<Vec<String>, (String, HashSet<String>)> = rest_primary_haplotypes.into_iter().take(hap_number-1).collect();
        let mut final_primary_haplotypes: HashMap<Vec<String>, (String, HashSet<String>)> = primary_haplotype_1.clone();
        for (path, (sequence, supported_reads)) in primary_haplotypes_2.iter() {
            final_primary_haplotypes.insert(path.clone(), (sequence.clone(), supported_reads.clone()));
        }
        final_primary_haplotypes
    } else {
        all_sequences
    };

    write_graph_path_fasta(&all_sequences_filtered, output_prefix);
    Ok(())
}