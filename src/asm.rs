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

#[derive(Debug, Clone)]
pub struct NodeInfo {
    pub seq: String,
    pub cigar: String,
    pub support_reads: String,
    pub allele_frequency: String,
    pub read_names: String,
}

#[derive(Debug, Clone)]
pub struct EdgeInfo {
    pub src: String,
    pub dst: String,
    pub overlap_ratio: f64,
    pub overlapping_reads: String,
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

pub fn find_germline_nodes(node_info: &HashMap<String, NodeInfo>, edge_info: &HashMap<String, Vec<String>>, hap_number: usize) -> (HashMap<String,NodeInfo>, HashMap<String, Vec<String>>) {
    let mut germline_nodes: HashMap<String, HashSet<String>> = HashMap::new();
    for (node, node_info) in node_info.iter() {
        let interval_name = node.split(".").collect::<Vec<_>>()[1];
        germline_nodes.entry(interval_name.to_string().clone()).or_insert(HashSet::new()).insert(node.clone());
    }

    let mut germline_nodes_sorted = HashMap::new();
    for (interval_name, node_vec) in germline_nodes.iter() {
        let mut node_support = Vec::new();
        for node in node_vec.iter() {
            let read_names = node_info.get(node).unwrap().read_names.clone();
            let read_names_list  = read_names.split(",").collect::<Vec<_>>();
            node_support.push((node.clone(), read_names_list.len()));
        }
        node_support.sort_by(|a, b| b.1.cmp(&a.1));

        // let noc_vec_len = &node_vec_.len();
        let n = node_support.len().min(hap_number);
        let node_vec_sorted = node_support.into_iter().take(n).collect::<Vec<_>>();
        for (node_id, _) in node_vec_sorted {
            germline_nodes_sorted.insert(node_id.clone(), node_info.get(&node_id).unwrap().clone());
        }    
    }
    let mut germline_edge_info = HashMap::new();
    for (g_node, g_node_info) in germline_nodes_sorted.iter() {
        if !edge_info.contains_key(g_node) {
            continue;
        }
        let next_nodes = edge_info.get(g_node).expect(&format!("Edge info not found for node: {}", g_node));
        for next_node in next_nodes {
            if germline_nodes_sorted.contains_key(next_node) {
                germline_edge_info.entry(g_node.clone()).or_insert(HashSet::new()).insert(next_node.clone());
            }
        }
    }
    let germline_edge_info_ = germline_edge_info.into_iter().map(|(k, v)| (k, v.into_iter().collect::<Vec<_>>())).collect();

    (germline_nodes_sorted, germline_edge_info_)
}

pub fn find_all_reads (node_info: &HashMap<String, NodeInfo>) -> HashSet<String> {
    let mut all_reads = HashSet::new();
    for (node, node_info) in node_info.iter() {
        let read_names = node_info.read_names.clone();
        all_reads.extend(read_names.split(",").map(|x| x.to_string()));
    }
    all_reads
}

/// Enumerate all paths in a directed acyclic graph using DFS
pub fn enumerate_all_paths(
    _node_info: &HashMap<String, NodeInfo>, 
    edge_info: &HashMap<String, Vec<String>>,
    read_constraint: bool
) -> std::result::Result<Vec<Vec<String>>, Box<dyn Error>> {
    let source_nodes = find_source_node(edge_info);
    let mut all_paths = Vec::new();
    let all_reads = find_all_reads(_node_info);
    if read_constraint {
        for src in source_nodes {
            let mut path = Vec::new();
            path.push(src.clone());
            dfs_traverse_with_read_constrains(&src, edge_info, _node_info, &mut path, &mut all_paths, &mut all_reads.clone());
        }

    }else{
        for src in source_nodes {
            let mut path = Vec::new();
            path.push(src.clone());
            dfs_traverse(&src, edge_info, &mut path, &mut all_paths);
        }

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
        // println!("current_node: {}, current_path: {:?}", current_node, current_path.len());
        all_paths.push(current_path.clone());
        return;
    }
    // println!("current_node: {}, current_path: {:?}", current_node, current_path.len();
    // Explore all outgoing edges
    if let Some(next_nodes) = edge_info.get(current_node) {
        for next_node in next_nodes {
            current_path.push(next_node.clone());
            dfs_traverse(next_node, edge_info, current_path, all_paths);
            current_path.pop(); // Backtrack
        }
    }
}


pub fn traverse_graph(node_info: &HashMap<String, NodeInfo>, edge_info: &HashMap<String, Vec<String>>, germline_only: bool, hap_number: usize, read_constraint: bool) -> std::result::Result<HashMap<Vec<String>, (String, HashSet<String>)>, Box<dyn Error>> {
    let all_paths = if germline_only {
        let (germline_nodes_sorted, germline_edge_info) = find_germline_nodes(node_info, edge_info, hap_number);
        enumerate_all_paths(&germline_nodes_sorted, &germline_edge_info, read_constraint)?
    } else {
        enumerate_all_paths(node_info, edge_info,read_constraint)?
    };
    // let all_paths = enumerate_all_paths(node_info, edge_info)?;
    
    // Get the sequence of each path
    let mut all_sequences = HashMap::new();
    for path in all_paths {
        let mut sequence = String::new();
        // let mut total_supported_reads = 0;
        let mut read_names = HashSet::new();
        for node in path.iter() {
            let node_info_dict = node_info.get(node).expect(&format!("Node {} not found in node_info, path Number: {:?}", node, path));
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
/// Recursive DFS to find all paths from a starting node with read constrains
fn dfs_traverse_with_read_constrains(
    current_node: &String,
    edge_info: &HashMap<String, Vec<String>>,
    node_info: &HashMap<String, NodeInfo>,
    current_path: &mut Vec<String>,
    all_paths: &mut Vec<Vec<String>>,
    read_intersection: &mut HashSet<String>,
) {
    // If no outgoing edges, this is a complete path
    if !edge_info.contains_key(current_node) {
        // println!("current_node: {}, current_path: {:?}", current_node, current_path.len());
        all_paths.push(current_path.clone());
        return;
    }
    if read_intersection.is_empty() {
        all_paths.push(current_path.clone());
        return;
    }
    // println!("current_node: {}, current_path: {:?}", current_node, current_path.len();
    // Explore all outgoing edges
    if let Some(next_nodes) = edge_info.get(current_node) {
        for next_node in next_nodes {
            current_path.push(next_node.clone());
            let read_names = node_info.get(next_node).unwrap().read_names.clone();
            let read_names_list = read_names.split(",").map(|x| x.to_string()).collect::<HashSet<_>>();
            let mut read_intersection_clone: HashSet<String> = read_intersection.intersection(&read_names_list).cloned().collect::<HashSet<_>>();
            dfs_traverse_with_read_constrains(next_node, edge_info, node_info, current_path, all_paths, &mut read_intersection_clone);
            current_path.pop(); // Backtrack
        }
    }
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
pub fn find_primary_haplotypes(all_sequences: &HashMap<Vec<String>, (String, HashSet<String>)>, haplotype_number: usize) -> HashMap<Vec<String>, (String, HashSet<String>)> {
    // sort all_sequences by the spanning length and the supported_reads
    let mut full_sequences = Vec::new();
    let mut all_reads = HashSet::new();
    for (path, (sequence, supported_reads)) in all_sequences.iter() {
        all_reads.extend(supported_reads.iter().cloned());
        let (start, end) = eval::find_alignment_intervals(path.iter().map(|x| x.as_str()).collect::<Vec<_>>()).unwrap();
        full_sequences.push((path.clone(), (sequence.clone(), supported_reads.clone(), end-start)));
    }
    full_sequences.sort_by(|a, b| b.1.1.len().cmp(&a.1.1.len()));


    // to be optimized
    let primary_haplotype_1: HashMap<Vec<String>, (String, HashSet<String>, usize)> = full_sequences.clone().into_iter().take(1).map(|(k, v)| (k.clone(), v.clone())).collect();
    let read_set_hap1 = primary_haplotype_1.values().map(|(_, supported_reads, _)| supported_reads).collect::<Vec<_>>()[0];
    let mut rest_primary_haplotypes: Vec<(Vec<String>, (String, HashSet<String>, usize))> = Vec::new();
    let read_set_hap2: HashSet<String> = (&all_reads - read_set_hap1).iter().cloned().collect();
    for (path, (sequence, supported_reads, span)) in full_sequences.clone().into_iter().skip(1) {
        let overlapped_reads: HashSet<String> = supported_reads.intersection(&read_set_hap2).cloned().collect();
        rest_primary_haplotypes.push((path.clone(), (sequence.clone(), overlapped_reads, span.clone())));
    }
    // sort the primary_haplotypes by the supported_reads
    rest_primary_haplotypes.sort_by(|a, b| b.1.1.len().cmp(&a.1.1.len()));
    let primary_haplotypes_2: HashMap<Vec<String>, (String, HashSet<String>, usize)> = rest_primary_haplotypes.into_iter().take(haplotype_number.min(full_sequences.len()) - 1).collect();
    let mut final_primary_haplotypes: HashMap<Vec<String>, (String, HashSet<String>)> = primary_haplotype_1.into_iter().map(|(k, v)| (k.clone(), (v.0.clone(), v.1.clone()))).collect();
    for (path, (sequence, supported_reads, span)) in primary_haplotypes_2.iter() {
        final_primary_haplotypes.insert(path.clone(), (sequence.clone(), supported_reads.clone()));
    }
    final_primary_haplotypes
}



pub fn start(graph_filename: &PathBuf, locus: &String, germline_only:bool, haplotype_number: usize, read_constraint: bool, output_prefix: &PathBuf) -> AnyhowResult<()> {
    let (node_info, edge_info) = load_graph(graph_filename).unwrap();
    // println!("node_info: {:?}", node_info);
    // println!("edge_info: {:?}", edge_info);
    let all_sequences = traverse_graph(&node_info, &edge_info, germline_only, haplotype_number, read_constraint).unwrap();
    info!("All sequences constructed: {}", all_sequences.len());
    let (ref_chromosome, ref_start, ref_end) = util::split_locus(locus.to_string());
    let mut hap_number = haplotype_number;
    let all_sequences_filtered = if germline_only {
        find_primary_haplotypes(&all_sequences, haplotype_number)
    } else {
        all_sequences
    };
    info!("All sequences filtered: {}", all_sequences_filtered.len());
    write_graph_path_fasta(&all_sequences_filtered, output_prefix);
    info!("All sequences written to fasta: {}", output_prefix.to_str().unwrap());
    Ok(())
}