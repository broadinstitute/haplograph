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
    pub support_reads: usize,
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
                support_reads: value["support_reads"].to_string().trim_matches('"').parse::<usize>().expect(&format!("Support reads not found for node: {}, {}", name, value["support_reads"].to_string())),
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

pub fn find_all_interval_names(node_info: &HashMap<String, NodeInfo>) -> HashMap<String, HashSet<String>> {
    let mut all_nodes: HashMap<String, HashSet<String>> = HashMap::new();
    for (node, node_info) in node_info.iter() {
        let interval_name = node.split(".").collect::<Vec<_>>()[1];
        all_nodes.entry(interval_name.to_string().clone()).or_insert(HashSet::new()).insert(node.clone());
    }
    all_nodes
}

pub fn find_germline_nodes(node_info: &HashMap<String, NodeInfo>, edge_info: &HashMap<String, Vec<String>>, hap_number: usize) -> (HashMap<String,NodeInfo>, HashMap<String, Vec<String>>, HashMap<String, HashSet<usize>>) {

    let all_nodes = find_all_interval_names(node_info);
    let mut all_nodes_sorted = all_nodes.into_iter().map(|(k, v)| (k, v.into_iter().collect::<Vec<_>>())).collect::<Vec<(String, Vec<String>)>>();
    all_nodes_sorted.sort_by(|a, b| a.0.cmp(&b.0));

    let mut germline_nodes_sorted = HashMap::new();
    let mut haplotype_reads:HashMap<usize, HashSet<String>> = HashMap::new();
    for i in 0.. hap_number {
        haplotype_reads.insert(i, HashSet::new());
    } // (hap_index, read_names) read specific to haplotype

    let mut node_haplotype: HashMap<String, HashSet<usize>> = HashMap::new(); // {node_names : hap} node specific to haplotype

    for (interval_name, node_vec) in all_nodes_sorted.iter() {
        let mut node_support = Vec::new();
        for node in node_vec.iter() {
            let read_name_list = get_read_name_list(node_info, node.clone());
            // let read_names = node_info.get(node).unwrap().read_names.clone();
            // let read_names_list  = read_names.split(",").collect::<Vec<_>>().iter().map(|x| x.to_string()).collect::<HashSet<_>>();
            node_support.push((node.clone(), read_name_list));
        }
        node_support.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

        // let noc_vec_len = &node_vec_.len();
        let n = node_support.len().min(hap_number);

        let node_vec_sorted = node_support.into_iter().take(n).collect::<Vec<_>>();
        for (hap_index, (node_id, read_names)) in node_vec_sorted.iter().enumerate() {
            germline_nodes_sorted.insert(node_id.clone(), node_info.get(&node_id.clone()).unwrap().clone());
            let mut unique_read_names = read_names.clone();
            // get unique read names from other parallele nodes
            for (hap_index_next, (node_id_next, read_names_next)) in node_vec_sorted.iter().enumerate() {
                if node_id_next == node_id {
                    continue;
                }
                unique_read_names = unique_read_names.difference(&read_names_next).cloned().collect::<HashSet<_>>();
            }

            // determine the haplotype index of reads
            if unique_read_names.is_empty(){
                continue;
            }  
            let mut max_overlap_reads = usize::MIN;
            let mut major_haplotype_index = hap_number + 1 as usize; 
            for (hap, read_names) in haplotype_reads.iter_mut() {
                let overlap_reads = read_names.intersection(&unique_read_names).cloned().collect::<HashSet<_>>();
                if overlap_reads.len() > max_overlap_reads {
                    max_overlap_reads = overlap_reads.len();
                    major_haplotype_index = hap.clone();
                }
            }

            if max_overlap_reads > 0 {
                node_haplotype.entry(node_id.clone()).or_insert(HashSet::new()).insert(major_haplotype_index.clone());
                haplotype_reads.entry(major_haplotype_index).or_insert(HashSet::new()).extend(unique_read_names.clone());
            }else{
                haplotype_reads.entry(hap_index).or_insert(HashSet::new()).extend(unique_read_names.clone());
                node_haplotype.entry(node_id.clone()).or_insert(HashSet::new()).insert(hap_index.clone());
            }
        }  
    }
    // println!("haplotype_reads: {:?}", haplotype_reads);

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

    (germline_nodes_sorted, germline_edge_info_, node_haplotype)
}

pub fn find_all_reads (node_info: &HashMap<String, NodeInfo>) -> HashSet<String> {
    let mut all_reads = HashSet::new();
    for (node, node_info) in node_info.iter() {
        let read_names = node_info.read_names.clone();
        all_reads.extend(read_names.split(",").map(|x| x.to_string()));
    }
    all_reads
}

pub fn get_read_name_list(node_info: &HashMap<String, NodeInfo>, node_id: String) -> HashSet<String> {
    let read_names = node_info.get(&node_id).unwrap().read_names.clone();
    let read_names_list = read_names.split(",").collect::<Vec<_>>().iter().map(|x| x.to_string()).collect::<HashSet<_>>();
    let read_names_list_clone = read_names_list.iter().map(|x| x.split("|").collect::<Vec<_>>()[0].to_string().clone()).collect::<HashSet<_>>();
    read_names_list_clone
}

pub fn find_best_node(node_info: &HashMap<String, NodeInfo>, current_nodes: &HashSet<String>) -> String {
    let mut best_node = String::new();
    let mut max_reads = 0;
    for node in current_nodes {
        let read_names = get_read_name_list(node_info, node.clone()); 
        if read_names.len() > max_reads {
            max_reads = read_names.len();
            best_node = node.clone();
        }
    }
    if best_node.is_empty() {
        info!("No best node found for haplotype: {:?}, {}", current_nodes, max_reads);
    }
    best_node
}

pub fn construct_sequences_from_haplotype(node_info: &HashMap<String, NodeInfo>, edge_info: &HashMap<String, Vec<String>>, node_haplotype: &HashMap<String, HashSet<usize>>, haplotype_number: usize, reference_sequence: String, reference_guided: bool) -> HashMap<usize, (Vec<String>, String, HashSet<String>)> {
    let mut haplotype_path = HashMap::new();
    let mut haplotype_sequences = HashMap::new();
    let mut haplotype_read_names = HashMap::new();
    let all_interval_nodes = find_all_interval_names(node_info);
    // find haplotype path
    for (node_id, haplotype_set) in node_haplotype.iter() {
        for haplotype in haplotype_set {
            haplotype_path.entry(haplotype.clone()).or_insert(Vec::new()).push(node_id.clone());
            let read_names_list_clone = get_read_name_list(node_info, node_id.clone());
            haplotype_read_names.insert(haplotype.clone(), read_names_list_clone.clone());
        }
    }

    // find haplotype sequence
    for (haplotype_num, node_list) in haplotype_path.iter() {
        let mut complete = true;
        let mut interval2nodes = HashMap::new();
        for node in node_list {
            let interval_name = node.split(".").collect::<Vec<_>>()[1];
            interval2nodes.entry(interval_name.to_string().clone()).or_insert(HashSet::new()).insert(node.clone());
        }
        let mut interval_list = interval2nodes.keys().collect::<Vec<_>>();
        interval_list.sort();
        let mut haplotype_sequence = String::new();
        let mut path: Vec<String> = Vec::new();
        let mut haplotype_reads_final = HashSet::new();
        for interval_name in interval_list.iter() {
            let current_nodes_by_assign = interval2nodes.get(*interval_name).unwrap();
            let current_nodes = if path.is_empty() {
                current_nodes_by_assign.clone()
            } else {
                let previous_node = path.last().unwrap().clone();
                let graph_current_nodes = if edge_info.contains_key(&previous_node) {
                    edge_info.get(&previous_node).unwrap().iter().map(|x| x.clone()).collect::<HashSet<_>>()
                } else {
                    current_nodes_by_assign.clone()
                };
                current_nodes_by_assign.intersection(&graph_current_nodes).cloned().collect::<HashSet<_>>()
            };
            
            if current_nodes.len() == 1 {
                let node = current_nodes.iter().next().unwrap();
                haplotype_sequence += &node_info.get(node).unwrap().seq.clone();
                path.push(node.clone());
                haplotype_reads_final.extend(get_read_name_list(node_info, node.clone()));
            } else if current_nodes.len() > 1 {
                let node = find_best_node(node_info, &current_nodes);
                haplotype_sequence += &node_info.get(&node).unwrap().seq.clone();
                path.push(node.clone());
                haplotype_reads_final.extend(get_read_name_list(node_info, node.clone()));
            }else{
                // no best node found, use the reference sequence
                info!("haplotype {} there is No connected node found for interval: {}, using reference sequence if guided by reference", haplotype_num, interval_name);
                if reference_guided{
                    let chromosome= interval_name.split(":").collect::<Vec<_>>()[0].to_string();
                    let position= interval_name.split(":").collect::<Vec<_>>()[1].split("-").collect::<Vec<_>>();
                    let start= position[0].parse::<usize>().unwrap() - 1;
                    let end= position[1].parse::<usize>().unwrap() - 1;
                    let reference_seq = reference_sequence[start..end].to_string();
                    haplotype_sequence += &reference_seq;
                }else{
                    complete = false;
                }
                
            }
        }
        // filter path based on edge_info
        if !reference_guided{
            for (i, node) in path.iter().enumerate() {
                if i == path.len() - 1 {
                    continue;
                }
                let next_node = path[i + 1].clone();
                if let Some(next_nodes) = edge_info.get(node) {
                    if !next_nodes.contains(&next_node) {
                        println!("haplotype_num: {}, Edge info not found for node: {}, next node: {}", haplotype_num, node, next_node);
                        complete = false;
                        break;
                    }
                }
            }
        }

        if complete{
            haplotype_sequences.insert(haplotype_num.clone(), (path.clone(), haplotype_sequence.clone(), haplotype_reads_final.clone()));
        }
    }

    haplotype_sequences
}

/// Enumerate all paths in a directed acyclic graph using DFS
pub fn enumerate_all_paths(
    _node_info: &HashMap<String, NodeInfo>, 
    edge_info: &HashMap<String, Vec<String>>,
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

pub fn enumerate_all_paths_with_haplotype(
    _node_info: &HashMap<String, NodeInfo>, 
    edge_info: &HashMap<String, Vec<String>>,
    node_haplotype: &HashMap<String, HashSet<usize>>,
    haplotype_number: usize,
) -> std::result::Result<Vec<Vec<String>>, Box<dyn Error>> {
    let source_nodes = find_source_node(edge_info);
    let mut all_paths = Vec::new();
    let mut haplotype_index = HashSet::new();
    for i in 0..haplotype_number {
        haplotype_index.insert(i);
    }
    for src in source_nodes {
        let mut path = Vec::new();
        path.push(src.clone());
        dfs_traverse_with_haplotype_constrains(&src, edge_info, &mut path, &mut all_paths, &mut haplotype_index, node_haplotype);
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
    // println!("current_node: {}, current_path: {:?}", current_node, current_path.len());
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

// /// Recursive DFS to find all paths from a starting node
fn dfs_traverse_with_haplotype_constrains(
    current_node: &String,
    edge_info: &HashMap<String, Vec<String>>,
    current_path: &mut Vec<String>,
    all_paths: &mut Vec<Vec<String>>,
    haplotype_index: &mut HashSet<usize>,
    node_haplotype: &HashMap<String, HashSet<usize>>
) {
    // println!("current_node: {}, read_intersection: {:?}", current_node, read_intersection.len());
    // If no outgoing edges, this is a complete path
    if !edge_info.contains_key(current_node) {
        if haplotype_index.len() > 0 {
            all_paths.push(current_path.clone());
        }
        return;
    }

    if haplotype_index.is_empty() {
        return;
    }

    let next_nodes = edge_info.get(current_node).unwrap();
    for next_node in next_nodes {
        let haplotypes_next = node_haplotype.get(&next_node.clone()).unwrap();
        let mut haplotype_intersection_clone = haplotype_index.intersection(&haplotypes_next).cloned().collect::<HashSet<_>>();
        current_path.push(next_node.clone());
        dfs_traverse_with_haplotype_constrains(next_node, edge_info,  current_path, all_paths, &mut haplotype_intersection_clone, node_haplotype);
        current_path.pop(); // Backtrack
    }
    
}

pub fn construct_sequences_from_haplotype_path(node_info: &HashMap<String, NodeInfo>, all_paths: &Vec<Vec<String>>) -> HashMap<usize, (Vec<String>, String, HashSet<String>)> {
    // Get the sequence of each path
    let mut all_sequences = HashMap::new();
    for (path_index, path) in all_paths.iter().enumerate() {
        let mut sequence = String::new();
        // let mut total_supported_reads = 0;
        let mut read_names = HashSet::new();
        for node in path.iter() {
            let node_info_dict = node_info.get(node).expect(&format!("Node {} not found in node_info, path Number: {:?}", node, path));
            let read_names_list_clone = get_read_name_list(node_info, node.clone());
            read_names.extend(read_names_list_clone.clone());
            let haplotype_seq = node_info_dict.seq.clone();
            sequence += &haplotype_seq;
            // total_supported_reads += supported_reads;
        }
        all_sequences.insert(path_index, (path.clone(), sequence, read_names.clone()));
    }
    all_sequences
}


pub fn write_graph_path_fasta(all_sequences: &HashMap<usize, (Vec<String>, String, HashSet<String>)>, output_filename: &PathBuf) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;
    let chars_per_line = 60;
    for (index,(path, sequence, supported_reads)) in all_sequences.iter(){
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
    // full_sequences.sort_by(|a, b| b.1.1.len().cmp(&a.1.1.len()));
    full_sequences.sort_by(|a, b| {
        // First compare reference spanning length
        b.1.2.cmp(&a.1.2)
            // Then compare supported_reads
            .then(b.1.1.len().cmp(&a.1.1.len()))
    });


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



pub fn start(graph_filename: &PathBuf, reference_fa: &String, germline_only:bool, haplotype_number: usize,  output_prefix: &PathBuf, reference_guided: bool) -> AnyhowResult<()> {
    let (node_info, edge_info) = load_graph(graph_filename).unwrap();
    // println!("node_info: {:?}", node_info.keys().collect::<Vec<_>>()[0]);
    info!("Reference guided: {}", reference_guided);
    let chromosome = node_info.keys().next().unwrap().split(".").collect::<Vec<_>>()[1].split(":").collect::<Vec<_>>()[0];
    let reference_chromosome = util::get_chromosome_ref_seq(&reference_fa, chromosome);
    let reference_sequence = String::from_utf8_lossy(reference_chromosome.iter().map(|x| x.seq()).collect::<Vec<_>>()[0]).to_string();
    

    info!("Traversing graph with germline_only: {}, hap_number: {}", germline_only, haplotype_number);
    let all_sequences = if germline_only {
        let (germline_nodes_sorted, germline_edge_info, node_haplotype) = find_germline_nodes(&node_info, &edge_info, haplotype_number);
        // enumerate_all_paths_with_haplotype(&germline_nodes_sorted, &germline_edge_info, &node_haplotype, hap_number)?
        construct_sequences_from_haplotype(&germline_nodes_sorted, &germline_edge_info, &node_haplotype, haplotype_number, reference_sequence, reference_guided)
    } else {
        let all_paths = enumerate_all_paths(&node_info, &edge_info).expect("Failed to enumerate all paths");
        construct_sequences_from_haplotype_path(&node_info, &all_paths)
    };

    info!("All sequences constructed: {}", all_sequences.len());

    let output_filename = PathBuf::from(format!("{}.fasta", output_prefix.to_string_lossy()));
    write_graph_path_fasta(&all_sequences, &output_filename);
    info!("All sequences written to fasta: {}", output_prefix.to_str().unwrap());
    Ok(())
}