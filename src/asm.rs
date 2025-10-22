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

pub fn find_parallele_nodes(node_info: &HashMap<String, NodeInfo>) -> HashMap<String, HashSet<String>> {
    let mut all_nodes: HashMap<String, HashSet<String>> = HashMap::new();
    for (node, node_info) in node_info.iter() {
        let interval_name = node.split(".").collect::<Vec<_>>()[1];
        all_nodes.entry(interval_name.to_string().clone()).or_insert(HashSet::new()).insert(node.clone());
    }
    all_nodes
}

pub fn identify_heterozygous_nodes(node_info: &HashMap<String, NodeInfo>, hap_number: usize, het_fold_threshold: f64) -> (HashMap<String, HashSet<String>>, HashMap<usize, HashSet<String>>) {
    let mut heterozygous_nodes: HashMap<String, HashSet<String>> = HashMap::new();
    let mut haplotype_reads: HashMap<usize, HashSet<String>> = HashMap::new();
    let all_nodes = find_parallele_nodes(node_info);
    for (interval_name, node_vec) in all_nodes.iter() {
        if node_vec.len() < 2{
            continue;
        }
        let mut node_support = Vec::new();
        for node in node_vec.iter() {
            let read_name_list = get_read_name_list(node_info, node.clone());
            node_support.push((node.clone(), read_name_list.len()));
        }
        node_support.sort_by(|a, b| b.1.cmp(&a.1));
        // println!("node_support: {:?}", node_support);
        let mut heterozygous = true;
        for i in 0..hap_number - 1 {
            let ratio = node_support[i].1 as f64 / node_support[i+1].1 as f64;
            // let total_reads = node_vec_sorted.iter().map(|(_, read_name_list_len)| read_name_list_len).sum::<usize>()/node_vec_sorted.len();
            if ratio > het_fold_threshold {
                heterozygous = false;
            }
        }
        if heterozygous {
            heterozygous_nodes.entry(interval_name.to_string().clone()).or_insert(HashSet::new()).extend(node_support.iter().take(hap_number).map(|(node, _)| node.clone()));
            // assign reads to haplotypes
            if haplotype_reads.is_empty() {
                for index in 0..hap_number {
                    let read_name_list = get_read_name_list(node_info, node_support[index].0.clone());
                    haplotype_reads.entry(index).or_insert(HashSet::new()).extend(read_name_list.clone());
                    // haplotype_reads[index] = haplotype_reads.get(index, set()) | find_read_name_list(node_info, node_support[index].0.clone());
                }
            }else{
                for index in 0..hap_number {
                    let read_name_list = get_read_name_list(node_info, node_support[index].0.clone());
                    let mut max_overlap_reads = 0;
                    let mut major_haplotype_index = hap_number + 1 as usize; 
                    for (hap, hap_reads) in haplotype_reads.iter_mut() {
                        let overlap_reads = hap_reads.intersection(&read_name_list).cloned().collect::<HashSet<_>>();
                        if overlap_reads.len() > max_overlap_reads {
                            max_overlap_reads = overlap_reads.len();
                            major_haplotype_index = hap.clone();
                        }
                    }
                    if max_overlap_reads > 0 {
                        haplotype_reads.entry(major_haplotype_index).or_insert(HashSet::new()).extend(read_name_list.clone());
                    }
                }
            }       
        }

    }
    (heterozygous_nodes, haplotype_reads)
}


pub fn assign_node_to_reads(node_info: &HashMap<String, NodeInfo>) -> HashMap<String, HashSet<String>> {
    let mut read_to_nodes = HashMap::new();
    let nodelist = node_info.keys().collect::<Vec<_>>();
    for node in nodelist {
        let read_names = get_read_name_list(node_info, node.clone());
        for read in read_names {
            read_to_nodes.entry(read).or_insert(HashSet::new()).insert(node.clone());
        }
    }
    read_to_nodes
}

pub fn assign_haplotype_nodes(node_info: &HashMap<String, NodeInfo>, haplotype_reads: &HashMap<usize, HashSet<String>>) -> HashMap<usize, HashSet<String>> {
    let read_to_nodes = assign_node_to_reads(node_info);
    let mut haplotype_nodes = HashMap::new();
    for (haplotype, reads) in haplotype_reads.iter() {
        for r in reads {
            haplotype_nodes.entry(haplotype.clone()).or_insert(HashSet::new()).extend(read_to_nodes.get(r).unwrap().iter().cloned());
        }
    }
    haplotype_nodes
}


pub fn find_unassigned_reads(node_info: &HashMap<String, NodeInfo>, haplotype_reads: &HashMap<usize, HashSet<String>>) -> HashSet<String> {
    let mut assigned_reads = HashSet::new();
    for reads in haplotype_reads.values() {
        assigned_reads.extend(reads.iter().cloned());
    }
    let all_reads = find_all_reads(node_info);
    all_reads.difference(&assigned_reads).cloned().collect()
}

pub fn assign_unassigned_reads(node_info: &HashMap<String, NodeInfo>, haplotype_reads: &HashMap<usize, HashSet<String>>) -> (HashMap<usize, HashSet<String>>, HashMap<usize, HashSet<String>>) {
    let read_to_nodes = assign_node_to_reads(node_info);
    let haplotype_nodes = assign_haplotype_nodes(node_info, haplotype_reads);
    let unassigned_reads = find_unassigned_reads(node_info, haplotype_reads);
    let mut haplotype_reads_new = haplotype_reads.clone();
    let mut haplotype_nodes_new = haplotype_nodes.clone();
    for read in unassigned_reads {
        let read_nodes = read_to_nodes.get(&read).unwrap();
        let mut max_overlap = 0;
        let mut max_haplotype = 100;
        for (haplotype, hap_nodes) in haplotype_nodes.iter() {
            let overlap_nodes = hap_nodes.intersection(read_nodes).cloned().collect::<HashSet<_>>();
            if overlap_nodes.len() > max_overlap {
                max_overlap = overlap_nodes.len();
                max_haplotype = haplotype.clone();
            }
        }
        if max_overlap > 0 {
            haplotype_reads_new.entry(max_haplotype).or_insert(HashSet::new()).insert(read.clone());
            haplotype_nodes_new.entry(max_haplotype).or_insert(HashSet::new()).extend(read_nodes.clone());
            
        }
    }
    (haplotype_reads_new, haplotype_nodes_new)
}

pub fn assign_haplotype_nodes_to_graph(node_info: &HashMap<String, NodeInfo>, haplotype_reads: &HashMap<usize, HashSet<String>>) -> HashMap<String, HashSet<usize>> {
    let mut node_haplotype = HashMap::new();
    let read_to_nodes = assign_node_to_reads(node_info);
    let all_nodes = find_parallele_nodes(node_info);
    if haplotype_reads.is_empty() {
        let haplotype = 0;
        for (interval_name, node_vec) in all_nodes.iter() {
            if node_vec.len() < 2{
                let node_id = node_vec.iter().next().unwrap().clone();
                node_haplotype.entry(node_id).or_insert(HashSet::new()).insert(haplotype);
            }else{
                let mut node_support = Vec::new();
                for node in node_vec.iter() {
                    let read_name_list = get_read_name_list(node_info, node.clone());
                    node_support.push((node.clone(), read_name_list.len()));
                }
                node_support.sort_by(|a, b| b.1.cmp(&a.1));
                let node_id = node_support[0].0.clone();
                node_haplotype.entry(node_id).or_insert(HashSet::new()).insert(haplotype);
            }
                
        }
    }else{
        for (haplotype, read_list) in haplotype_reads.iter() {
            for read in read_list {
                let node_list = read_to_nodes.get(read).unwrap();
                for node in node_list {
                    node_haplotype.entry(node.clone()).or_insert(HashSet::new()).insert(haplotype.clone());
                }
            }
        }
    }

    node_haplotype
}

pub fn assign_haplotype_to_nodes(haplotype_nodes: &HashMap<usize, HashSet<String>>) -> HashMap<String, HashSet<usize>> {
    let mut node_haplotype = HashMap::new();
    for (haplotype, nodes) in haplotype_nodes.iter() {
        for node in nodes {
            node_haplotype.entry(node.clone()).or_insert(HashSet::new()).insert(haplotype.clone());
        }
    }
    node_haplotype
}




pub fn find_all_reads (node_info: &HashMap<String, NodeInfo>) -> HashSet<String> {
    let mut all_reads = HashSet::new();
    for (node, node_infomation) in node_info.iter() {
        let read_names = get_read_name_list(&node_info, node.clone());
        all_reads.extend(read_names.iter().cloned());
    }
    all_reads
}

pub fn get_read_name_list(node_info: &HashMap<String, NodeInfo>, node_id: String) -> HashSet<String> {
    let read_names = node_info.get(&node_id).unwrap().read_names.clone();
    let read_names_list = read_names.split(",").collect::<Vec<_>>().iter().map(|x| x.to_string()).collect::<HashSet<_>>();
    let read_names_list_clone = read_names_list.iter().map(|x| x.split("|").collect::<Vec<_>>()[0].to_string().clone()).collect::<HashSet<_>>();
    read_names_list_clone
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

// /// Recursive DFS to find all paths from a starting node
fn dfs_traverse_with_haplotype_constrains(
    current_node: &String,
    edge_info: &HashMap<String, Vec<String>>,
    current_path: &mut Vec<String>,
    all_paths: &mut Vec<Vec<String>>,
    haplotype_index: &mut HashSet<usize>,
    node_haplotype: &HashMap<String, HashSet<usize>>
) {

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
        let haplotypes_next= if (node_haplotype.contains_key(next_node)) {
            node_haplotype.get(&next_node.clone()).unwrap()
        } else {
            &HashSet::new()
        };
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


pub fn write_graph_path_fasta(all_sequences: &Vec<(Vec<String>, String, HashSet<String>, usize)>, output_filename: &PathBuf) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;
    let chars_per_line = 60;
    for (index,(path, sequence, supported_reads, span)) in all_sequences.iter().enumerate(){
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
pub fn find_full_range_haplotypes(node_info: &HashMap<String, NodeInfo>, all_sequences: &HashMap<usize, (Vec<String>, String, HashSet<String>)>) -> Vec<(Vec<String>, String, HashSet<String>, usize)> {
    // sort all_sequences by the spanning length and the supported_reads
    let node_list = node_info.keys().collect::<Vec<_>>();
    let full_range = eval::find_alignment_intervals(node_list.iter().map(|x| x.as_str()).collect::<Vec<_>>()).unwrap();
    let mut full_sequences = Vec::new();
    let mut all_reads = HashSet::new();
    for (hap_index,(path, sequence, supported_reads)) in all_sequences.iter() {
        all_reads.extend(supported_reads.iter().cloned());
        let (start, end) = eval::find_alignment_intervals(path.iter().map(|x| x.as_str()).collect::<Vec<_>>()).unwrap();
        if start <= full_range.0 && end >= full_range.1 {
            full_sequences.push((path.clone(), sequence.clone(), supported_reads.clone(), end-start));
        }
    }
    full_sequences.sort_by(|a, b| b.2.len().cmp(&a.2.len()));
    full_sequences
}

pub fn find_parallele_nodes_from_nodelist(nodelist: &Vec<String>) -> HashMap<String, HashSet<String>> {
    let mut interval_node = HashMap::new();
    for node in nodelist {
        let interval = node.split(".").collect::<Vec<_>>()[1];
        interval_node.entry(interval.to_string()).or_insert(HashSet::new()).insert(node.clone());
    }
    interval_node
}

pub fn filter_haplotype_nodes(node_info: &HashMap<String, NodeInfo>, haplotype_nodes: &HashMap<usize, HashSet<String>>, haplotype_reads: &HashMap<usize, HashSet<String>>) -> HashMap<usize, HashSet<String>> {
    let mut filtered_haplotype_nodes = HashMap::new();
    for (hap, nodes_list) in haplotype_nodes.iter() {
        let interval_node = find_parallele_nodes_from_nodelist(&nodes_list.iter().cloned().collect::<Vec<_>>());
        let hap_specific_reads = haplotype_reads.get(hap).unwrap();
        for (interval, nodes) in interval_node.iter() {
            if nodes.len() > 1{
                let mut best_node = "".to_string();
                let mut best_read_count = 0;
                for node in nodes {
                    let read_names = get_read_name_list(node_info, node.clone()).intersection(&hap_specific_reads.clone()).cloned().collect::<HashSet<_>>();
                    if read_names.len() > best_read_count {
                        best_node = node.clone();
                        best_read_count = read_names.len();
                    }
                }
                if best_read_count > 0 {
                    filtered_haplotype_nodes.entry(hap.clone()).or_insert(HashSet::new()).insert(best_node.clone());
                }
            }else if nodes.len() == 1{
                filtered_haplotype_nodes.entry(hap.clone()).or_insert(HashSet::new()).insert(nodes.iter().next().unwrap().clone());
            }else{
                warn!("hap: {}, interval: {}, nodes: {:?}", hap, interval, nodes);
            }   

        }
    }
    filtered_haplotype_nodes
}

pub fn find_node_haplotype(node_info: &HashMap<String, NodeInfo>, hap_number: usize, het_fold_threshold: f64) -> (HashMap<usize, HashSet<String>>, HashMap<String, HashSet<usize>>) {
    let (heterozygous_nodes, haplotype_reads) = identify_heterozygous_nodes(node_info, hap_number, het_fold_threshold);
    println!("heterozygous_nodes: {:?}", heterozygous_nodes);
    println!("haplotype_reads: {:?}", haplotype_reads);
    let haplotype_nodes = assign_haplotype_nodes(node_info, &haplotype_reads);
    let read_to_nodes = assign_node_to_reads(node_info);
    let (haplotype_reads_new, haplotype_nodes_new) = assign_unassigned_reads(node_info, &haplotype_reads);
    // let node_haplotype = assign_haplotype_nodes_to_graph(node_info, &haplotype_reads_new);
    if haplotype_reads_new.is_empty() {
        let node_haplotype = assign_haplotype_nodes_to_graph(node_info, &haplotype_reads_new);
        return (haplotype_reads_new, node_haplotype);
    }else{
        let filtered_haplotype_nodes = filter_haplotype_nodes(node_info, &haplotype_nodes_new, &haplotype_reads_new);
        let node_haplotype = assign_haplotype_to_nodes(&filtered_haplotype_nodes);
        return (haplotype_reads_new, node_haplotype);
    }
}

pub fn start(graph_filename: &PathBuf, germline_only:bool, haplotype_number: usize,  output_prefix: &PathBuf, het_fold_threshold: f64) -> AnyhowResult<()> {
    let (node_info, edge_info) = load_graph(graph_filename).unwrap();
    info!("Traversing graph with germline_only: {}, hap_number: {}", germline_only, haplotype_number);

    let (haplotype_reads, node_haplotype) = find_node_haplotype(&node_info, haplotype_number, het_fold_threshold);
    let all_paths = enumerate_all_paths_with_haplotype(&node_info, &edge_info, &node_haplotype, haplotype_number).expect("Failed to enumerate all paths");
    let allseq = construct_sequences_from_haplotype_path(&node_info, &all_paths);
    let primary_haplotypes = find_full_range_haplotypes(&node_info, &allseq);
    info!("All sequences constructed: {}", primary_haplotypes.len());
    let output_filename = PathBuf::from(format!("{}.fasta", output_prefix.to_string_lossy()));
    let _ =write_graph_path_fasta(&primary_haplotypes, &output_filename);
    info!("All sequences written to fasta: {}", output_prefix.to_str().unwrap());
    Ok(())
}