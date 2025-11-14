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
use std::path::Path;
use ndarray::{Array2, Array1};

#[derive(Debug, Clone)]
pub struct NodeInfo {
    pub seq: String,
    pub cigar: String,
    pub support_reads: usize,
    pub allele_frequency: String,
    pub read_names: String,
    pub methyl_info: HashMap<usize, f32>,
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
            let methyl_info = value["mod_score_dict"].as_str().unwrap_or_default();
            let methyl_info_list = if methyl_info.is_empty() { Vec::new() } else { methyl_info.split(",").collect::<Vec<_>>() };
            let methyl_info_dict = if methyl_info_list.len() > 0 { methyl_info_list.iter().map(|x| x.split(":").collect::<Vec<_>>()).collect::<Vec<_>>().iter().map(|x| (x[0].parse::<usize>().unwrap(), x[1].parse::<f32>().unwrap())).collect::<HashMap<usize, f32>>() } else { HashMap::new() };
            // let methyl_info_dict = methyl_info.iter().map(|x| x.split(":").collect::<Vec<_>>()).collect::<Vec<_>>().iter().map(|x| (x[0].parse::<usize>().unwrap(), x[1].parse::<f32>().unwrap())).collect::<HashMap<usize, f32>>();
            let value_info = NodeInfo {
                seq: seq.to_string(),
                cigar: value["cigar"].to_string(),
                support_reads: value["support_reads"].to_string().trim_matches('"').parse::<usize>().expect(&format!("Support reads not found for node: {}, {}", name, value["support_reads"].to_string())),
                allele_frequency: value["allele_frequency"].to_string(),
                read_names: value["read_names"].to_string(),
                methyl_info: methyl_info_dict,
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

pub fn identify_heterozygous_nodes(node_info: &HashMap<String, NodeInfo>, hap_number: usize, het_fold_threshold: f64) -> HashMap<String, HashSet<String>> {
    let mut heterozygous_nodes: HashMap<String, HashSet<String>> = HashMap::new();
    let all_nodes = find_parallele_nodes(node_info);
    let mut interval_list = all_nodes.keys().collect::<Vec<_>>();
    interval_list.sort_by(|a, b| util::split_locus(a.to_string()).1.cmp(&util::split_locus(b.to_string()).1)); // ascending order a < b
    for interval_name in interval_list.iter() {
        let node_vec = all_nodes.get(interval_name.clone()).unwrap();
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
            if ratio > het_fold_threshold {
                heterozygous = false;
                break;
            }
        }
        if heterozygous {
            heterozygous_nodes.entry(interval_name.to_string().clone()).or_insert(HashSet::new()).extend(node_support.iter().take(hap_number).map(|(node, _)| node.clone()));
        }
    }
    heterozygous_nodes
}

pub fn assign_haplotype_reads(node_info: &HashMap<String, NodeInfo>, heterozygous_nodes: &HashMap<String, HashSet<String>>, hap_number: usize) -> HashMap<usize, HashSet<String>> {
    // assign reads to haplotypes
    let mut haplotype_reads = HashMap::new();
    let mut interval_list = heterozygous_nodes.keys().collect::<Vec<_>>();
    interval_list.sort_by(|a, b|util::split_locus(a.to_string()).1.cmp(&util::split_locus(b.to_string()).1));
    for interval_name in interval_list.iter() {
        let mut node_set = heterozygous_nodes.get(interval_name.clone()).unwrap().clone().iter().cloned().collect::<Vec<_>>();
        node_set.sort_by(|a, b| get_read_name_list(node_info, a.clone()).len().cmp(&get_read_name_list(node_info, b.clone()).len()));
        if node_set.len() < hap_number{
            continue;
        }
        if haplotype_reads.is_empty() {
            for (index, node) in node_set.iter().enumerate() {
                println!("node: {:?}, index: {:?}", node, index);
                let read_name_list = get_read_name_list(node_info, node.clone());
                haplotype_reads.entry(index).or_insert(HashSet::new()).extend(read_name_list.clone());
             }
        }else{
            for (index, node) in node_set.iter().enumerate() {
                let read_name_list = get_read_name_list(node_info, node.clone());
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
    
    haplotype_reads
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

pub fn construct_heterozygous_nodes_matrix(node_info: &HashMap<String, NodeInfo>, node_list:Vec<String>) -> Array2<f64> {
    let mut read_vec = HashSet::new();
    for node in node_list.iter() {
        let read_name_list = get_read_name_list(node_info, node.clone());
        read_vec.extend(read_name_list.iter().cloned());
    } 
    let read_list = read_vec.iter().cloned().collect::<Vec<_>>();
    let mut matrix = Array2::<f64>::zeros((node_list.len(), read_list.len()));
    for (n_index, node_name) in node_list.iter().enumerate() {
        let read_name_list = get_read_name_list(node_info, node_name.clone());
        for read in read_name_list {
            let read_index = read_list.iter().position(|x| x == &read).unwrap();
            matrix[[n_index, read_index]] = 1.0;
        }
    }
    matrix 
}

pub fn filter_heterozygous_nodes(node_info: &HashMap<String, NodeInfo>, heterozygous_nodes: &HashMap<String, HashSet<String>>) -> HashMap<String, HashSet<String>> {
    let mut nodelist = HashSet::new();
    for (interval_name, node_set) in heterozygous_nodes.iter() {
        nodelist.extend(node_set.iter().cloned());
    }
    let node_list = nodelist.iter().cloned().collect::<Vec<_>>();
    let matrix = construct_heterozygous_nodes_matrix(node_info, node_list.clone());
    let node_list_filtered = util::permutation_test(&matrix, 0.1, 100, node_list.clone());
    let mut filtered_heterozygous_nodes = HashMap::new();
    for node_id in node_list_filtered.iter() {
        let interval_name = node_id.split(".").collect::<Vec<_>>()[1].to_string();
        filtered_heterozygous_nodes.entry(interval_name.clone()).or_insert(HashSet::new()).insert(node_id.clone());
    }
    filtered_heterozygous_nodes
}

pub fn assign_unassigned_reads(node_info: &HashMap<String, NodeInfo>, haplotype_reads: &HashMap<usize, HashSet<String>>) -> (HashMap<usize, HashSet<String>>, HashMap<usize, HashSet<String>>) {
    let haplotype_nodes = assign_haplotype_nodes(node_info, haplotype_reads);
    let read_to_nodes = assign_node_to_reads(node_info);
    let unassigned_reads = find_unassigned_reads(node_info, haplotype_reads);
    let mut haplotype_reads_new = haplotype_reads.clone();
    let mut haplotype_nodes_new = haplotype_nodes.clone();
    for read in unassigned_reads {
        let read_nodes = read_to_nodes.get(&read).unwrap();
        // let mut max_overlap = 0;
        // let mut max_haplotype = 100;
        let mut find_haplotype = false;
        for (haplotype, hap_nodes) in haplotype_nodes.iter() {
            let overlap_nodes = hap_nodes.intersection(read_nodes).cloned().collect::<HashSet<_>>();
            if overlap_nodes.len() > 0 {
                haplotype_reads_new.entry(haplotype.clone()).or_insert(HashSet::new()).insert(read.clone());
                haplotype_nodes_new.entry(haplotype.clone()).or_insert(HashSet::new()).extend(read_nodes.clone());
                find_haplotype = true;
            }
        }
        // homologous reads
        if !find_haplotype {
            for hap in haplotype_reads_new.clone().keys() {
                haplotype_reads_new.entry(hap.clone()).or_insert(HashSet::new()).insert(read.clone());
                haplotype_nodes_new.entry(hap.clone()).or_insert(HashSet::new()).extend(read_nodes.clone());
            }
        }
        // if max_overlap > 0 {
        //     haplotype_reads_new.entry(max_haplotype).or_insert(HashSet::new()).insert(read.clone());
        //     haplotype_nodes_new.entry(max_haplotype).or_insert(HashSet::new()).extend(read_nodes.clone());
        // }else{
        //     for hap in haplotype_reads_new.clone().keys() {
        //         haplotype_reads_new.entry(*hap).or_insert(HashSet::new()).insert(read.clone());
        //         haplotype_nodes_new.entry(*hap).or_insert(HashSet::new()).extend(read_nodes.clone());
        //     }
        // }
    }
    (haplotype_reads_new, haplotype_nodes_new)
}



pub fn find_most_supported_path(node_info: &HashMap<String, NodeInfo>, edge_info: &HashMap<String, Vec<String>>) -> HashMap<String, HashSet<usize>> {
    let mut node_haplotype = HashMap::new();
    let all_nodes = find_parallele_nodes(node_info);
    let mut interval_list = all_nodes.keys().collect::<Vec<_>>();
    interval_list.sort_by(|a, b| a.cmp(b)); // ascending order a < b
    // should also consider the edge connectivity
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
            node_support.sort_by(|a, b| b.1.cmp(&a.1)); // descending order a > b
            let node_id = node_support[0].0.clone();
            node_haplotype.entry(node_id).or_insert(HashSet::new()).insert(haplotype);
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
    let read_names_list_clone = read_names_list.iter().map(|x| x.split("|").collect::<Vec<_>>()[0].to_string().replace("\"", "").clone()).collect::<HashSet<_>>();
    read_names_list_clone
}

pub fn enumerate_all_paths_with_haplotype(
    _node_info: &HashMap<String, NodeInfo>, 
    edge_info: &HashMap<String, Vec<String>>,
    node_haplotype: &HashMap<String, HashSet<usize>>,
    haplotype_number: usize,
) -> std::result::Result<Vec<(Vec<String>, HashSet<usize>)>, Box<dyn Error>> {
    let source_nodes = find_source_node(edge_info);
    let mut all_paths = Vec::new();
    let mut haplotype_index = HashSet::new();
    for i in 0..haplotype_number {
        haplotype_index.insert(i);
    }

    for src in source_nodes.clone() {
        let mut haplotype_intersection_src = node_haplotype.get(&src).unwrap_or(&HashSet::new()).intersection(&haplotype_index).cloned().collect::<HashSet<_>>();
        if haplotype_intersection_src.len() > 0 {
            let mut path = Vec::new();
            path.push(src.clone());
            dfs_traverse_with_haplotype_constrains(&src, edge_info, &mut path, &mut all_paths, &mut haplotype_intersection_src, node_haplotype);
        }
    }

    
    Ok(all_paths)
}

// /// Recursive DFS to find all paths from a starting node
fn dfs_traverse_with_haplotype_constrains(
    current_node: &String,
    edge_info: &HashMap<String, Vec<String>>,
    current_path: &mut Vec<String>,
    all_paths: &mut Vec<(Vec<String>, HashSet<usize>)>,
    haplotype_index: &mut HashSet<usize>,
    node_haplotype: &HashMap<String, HashSet<usize>>
) {

    if !edge_info.contains_key(current_node) {
        if haplotype_index.len() > 0 {
            all_paths.push((current_path.clone(), haplotype_index.clone()));
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

pub fn construct_sequences_from_haplotype_path(node_info: &HashMap<String, NodeInfo>, all_paths: &Vec<(Vec<String>, HashSet<usize>)>) -> HashMap<usize, Vec<(Vec<String>, String, HashSet<String>)>> {
    // Get the sequence of each path
    let mut all_sequences = HashMap::new();
    for (path_index, (path, haplotype_index)) in all_paths.iter().enumerate() {
        let mut sequence = String::new();
        let mut read_names = HashSet::new();
        
        for node in path.iter() {
            let node_info_dict = node_info.get(node).expect(&format!("Node {} not found in node_info, path Number: {:?}", node, path));
            let read_names_list_clone = get_read_name_list(node_info, node.clone());
            read_names.extend(read_names_list_clone.clone());
            let haplotype_seq = node_info_dict.seq.clone();
            sequence += &haplotype_seq;
            // total_supported_reads += supported_reads;
        }
        for hap_ind in haplotype_index.iter() {
            info!("path: {:?}, haplotype_index: {:?}, sequence length: {:?}", path.len(), hap_ind, sequence.len());
            all_sequences.entry(hap_ind.clone()).or_insert(Vec::new()).push((path.clone(), sequence.clone(), read_names.clone()));
        }
    }
    all_sequences
}


pub fn write_graph_path_fasta(all_sequences: &HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)>, output_filename: &PathBuf) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;
    let chars_per_line = 60;
    for (index,(path, sequence, supported_reads, supports, span)) in all_sequences.iter(){
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

pub fn get_supports (node_info: &HashMap<String, NodeInfo>, path: &Vec<String>) -> usize {
    let mut supports = 0 as usize;
    for node in path.iter() {
        let read_set = get_read_name_list(node_info, node.clone());
        supports += read_set.len();
    }
    supports
}

pub fn find_full_range_haplotypes(node_info: &HashMap<String, NodeInfo>,node_haplotype: &HashMap<String, HashSet<usize>>, all_sequences: &HashMap<usize, Vec<(Vec<String>, String, HashSet<String> )>>, hap_number: usize) -> HashMap<usize, (Vec<String>, String, HashSet<String>,usize, usize)> {
    // sort all_sequences by the spanning length and the supported_reads
    let node_list = node_info.keys().collect::<Vec<_>>();
    // let full_range = eval::find_alignment_intervals(node_list.iter().map(|x| x.as_str()).collect::<Vec<_>>()).unwrap();
    let mut best_paths = HashMap::new();
    for (hap_index,path_list) in all_sequences.iter() {
        // for each haplotype, select the best path
        let mut full_sequences = Vec::new();
        for (index, (path, sequence, supported_reads)) in path_list.iter().enumerate() {
            let (start, end) = eval::find_alignment_intervals(path.iter().map(|x| x.as_str()).collect::<Vec<_>>()).unwrap();
            let supports = get_supports(node_info, path);
            info!("path: {:?}, supports: {:?}, supported_reads: {:?}, span: {:?}", hap_index, supports, supported_reads.len(), end-start );

            full_sequences.push((path.clone(), sequence.clone(), supported_reads.clone(), supports, end-start));
            
        }
        //first compare the 4th element, then the 3th element
        // full_sequences.sort_by(|a, b| b.3.cmp(&a.3).then(b.4.cmp(&a.4)));
        full_sequences.sort_by(|a, b| b.2.len().cmp(&a.2.len()).then(b.3.cmp(&a.3)));
        best_paths.insert(hap_index.clone(),full_sequences[0].clone());
    }
    best_paths
    
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
    
    for (hap, node_list) in haplotype_nodes.iter() {
        let interval_node = find_parallele_nodes_from_nodelist(&node_list.iter().map(|x| x.clone()).collect::<Vec<_>>());
        let mut interval_list = interval_node.keys().collect::<Vec<_>>();
        interval_list.sort_by(|a, b| a.cmp(b)); // ascending order a < b
    
        let read_list = haplotype_reads.get(hap).unwrap().clone();
        for interval in interval_list.iter() {
            let nodes = interval_node.get(*interval).unwrap().clone();
            let mut best_node = "".to_string();
            let mut best_read_count = 0;
            for node in nodes.clone() {
                let node_read_names = get_read_name_list(node_info, node.clone());
                let intersection_count = read_list.intersection(&node_read_names).count();
                if intersection_count > best_read_count {
                    best_node = node.clone();
                    best_read_count = intersection_count;
                }
            }
            if best_read_count > 0 {
                filtered_haplotype_nodes.entry(hap.clone()).or_insert(HashSet::new()).insert(best_node.clone());
            }else{
                warn!("hap: {}, interval: {}, nodes: {:?}", hap, interval, nodes.clone().iter().map(|x| x.clone()).collect::<Vec<_>>().join(", "));
                // let mut best_node = "".to_string();
                // let mut best_read_count = 0;
                // for node in nodes.clone() {
                //     let node_read_names = get_read_name_list(node_info, node.clone());
                //     if node_read_names.len() > best_read_count {
                //         best_node = node.clone();
                //         best_read_count = node_read_names.len();
                //     }
                // }
                // filtered_haplotype_nodes.entry(hap.clone()).or_insert(HashSet::new()).insert(best_node.clone());
            }
        }
    }
    filtered_haplotype_nodes
}

pub fn find_node_haplotype(node_info: &HashMap<String, NodeInfo>, edge_info: &HashMap<String, Vec<String>>, hap_number: usize, het_fold_threshold: f64) -> (HashMap<usize, HashSet<String>>, HashMap<String, HashSet<usize>>) {
    let heterozygous_nodes = identify_heterozygous_nodes(node_info, hap_number, het_fold_threshold);
    info!("heterozygous_nodes: {:?}", heterozygous_nodes.len());
    let filtered_heterozygous_nodes = filter_heterozygous_nodes(node_info, &heterozygous_nodes);
    info!("filtered_heterozygous_nodes: {:?}", filtered_heterozygous_nodes.len());
    
    let haplotype_reads = assign_haplotype_reads(node_info, &filtered_heterozygous_nodes, hap_number);
    
    let (haplotype_reads_new, haplotype_nodes_new) = assign_unassigned_reads(node_info, &haplotype_reads);
    info!("haplotype_reads: {:?}", haplotype_reads.iter().map(|(hap, reads)| format!("hap: {}, reads: {}", hap, reads.len())).collect::<Vec<_>>().join(", "));
    info!("haplotype_reads_new: {:?}", haplotype_reads_new.iter().map(|(hap, reads)| format!("hap: {}, reads: {}", hap, reads.len())).collect::<Vec<_>>().join(", "));
    
    if haplotype_reads_new.is_empty() {
        let node_haplotype = find_most_supported_path(node_info, &edge_info);
        return (haplotype_reads_new, node_haplotype);
    }else{
        let filtered_haplotype_nodes = filter_haplotype_nodes(node_info, &haplotype_nodes_new, &haplotype_reads_new);
        let node_haplotype = assign_haplotype_to_nodes(&filtered_haplotype_nodes);
        // println!("node_haplotype: {:?}", node_haplotype);
        return (haplotype_reads_new, node_haplotype);
    }
}

/// Maps positions from alternate sequence to reference sequence based on CIGAR string
pub fn mapping_to_reference_coordinates(cigar: &str, ref_start: usize) -> HashMap<usize, usize> {
    let mut ref_pos = 0;
    let mut alt_pos = 0;
    
    // Parse CIGAR string into operations
    let mut operations = Vec::new();
    let mut num = String::new();
    
    for c in cigar.chars() {
        if c.is_digit(10) {
            num.push(c);
        } else {
            if !num.is_empty() {
                let length = num.parse::<usize>().unwrap();
                operations.push((length, c));
                num.clear();
            }
        }
    }
    
    // Process each operation
    let mut position_mapping = HashMap::new();
    
    for (length, op) in operations {
        match op {
            '=' | 'M' => {  // Match
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    position_mapping.insert(alt_pos + i, pos);
                }
                ref_pos += length;
                alt_pos += length;
            },
            'X' => {  // Mismatch
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    position_mapping.insert(alt_pos + i, pos);
                }
                ref_pos += length;
                alt_pos += length;
            },
            'I' => {  // Insertion
                let pos = ref_start + ref_pos - 1;  // Mapping to one base pair before
                for i in 0..length {
                    position_mapping.insert(alt_pos + i, pos);
                }
                alt_pos += length;
            },
            'D' => {  // Deletion, nothing will map back to reference coordinates
                ref_pos += length;
            },
            _ => {
                warn!("Unexpected CIGAR operation: {}", op);
            }
        }
    }
    
    position_mapping
}

fn write_methyl_bed(
    methyl_info: &HashMap<(usize, usize), (f32, usize)>,
    output_prefix: &PathBuf,
    haplotype_index: usize,
    chromosome: &str,
) -> std::io::Result<()> {
    let output_file = output_prefix.with_extension(format!("Hap.{}.bed", haplotype_index));
    let mut file = File::create(Path::new(&output_file))?;

    // Write BED header
    writeln!(file, "##fileformat=BED")?;
    writeln!(file, "##haplotype={}", haplotype_index + 1)?;
    writeln!(
        file,
        "#CHROM\tRef_start\tRef_end\tMod_rate\tAsm_start\tAsm_end\tMotif\tCoverage"
    )?;
    // let min_prob = 0.5;
    // let mut methylation_signal: HashMap<(usize, Option<usize>), (f64, f64)> = HashMap::new();
    let mut methyl_info_vec = Vec::new();
    for ((ref_pos, asm_pos), (score, coverage)) in methyl_info.iter() {
        // 1-based coordinates
        let ref_pos_start = ref_pos + 1;
        let ref_pos_end = ref_pos + 2;
        let asm_pos_start = asm_pos + 1;
        let asm_pos_end = asm_pos + 2;
        let motif = "CG".to_string();
        methyl_info_vec.push((chromosome, ref_pos_start, ref_pos_end, asm_pos_start, asm_pos_end, motif, score, coverage));

    }
    methyl_info_vec.sort_by_key(|item| item.1);
    for methyl_list in methyl_info_vec.iter(){
        let (chromosome, ref_pos_start, ref_pos_end, asm_pos_start, asm_pos_end, motif, methyl_rate, coverage) = methyl_list;
        writeln!(file, 
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chromosome, ref_pos_start, ref_pos_end,  methyl_rate, asm_pos_start, asm_pos_end, motif, coverage
            )?;

    }

    Ok(())
}


pub fn call_methylation(node_info: &HashMap<String, NodeInfo>, all_paths: HashMap<usize, (Vec<String>, String, HashSet<String>,usize, usize)>, output_prefix: &PathBuf) -> () {
    // call methylation
    for (hap_index, (path, sequence, supported_reads, supports, span)) in all_paths.iter() {
        let mut methyl_info_dict = HashMap::new();
        let mut spos = 0;
        let chromosome = path[0].split(".").collect::<Vec<_>>()[1].split(":").collect::<Vec<_>>()[0];
        for node in path.iter() {
            let methyl_info = node_info.get(node).unwrap().methyl_info.clone();
            let cigar = node_info.get(node).unwrap().cigar.clone();
            let ref_start = eval::find_alignment_intervals(vec![node.clone()].iter().map(|x| x.as_str()).collect::<Vec<_>>()).unwrap().0;
            let position_mapping = mapping_to_reference_coordinates(&cigar, ref_start);
            // println!("position_mapping: {:?}", position_mapping);
            let node_seq = node_info.get(node).unwrap().seq.clone();
            let node_coverage = node_info.get(node).unwrap().support_reads;
            for (pos, score) in methyl_info.into_iter() {
                let asm_pos = pos + spos;
                let ref_pos  = position_mapping.get(&pos).unwrap_or(&0).clone();
                // println!("ref_pos: {}, asm_pos: {}, score: {}", ref_pos, asm_pos, score);
                methyl_info_dict.insert((ref_pos, asm_pos), (score, node_coverage));
            }
            spos += node_seq.len();
        }
        let _ = write_methyl_bed(&methyl_info_dict, &output_prefix, *hap_index,  chromosome);
        
    }
}

pub fn start(graph_filename: &PathBuf, germline_only:bool, haplotype_number: usize,  output_prefix: &PathBuf, het_fold_threshold: f64) -> AnyhowResult<()> {
    let (node_info, edge_info) = load_graph(graph_filename).unwrap();
    info!("Traversing graph with germline_only: {}, hap_number: {}", germline_only, haplotype_number);

    let (haplotype_reads, node_haplotype) = find_node_haplotype(&node_info, &edge_info, haplotype_number, het_fold_threshold);
    let all_paths = enumerate_all_paths_with_haplotype(&node_info, &edge_info, &node_haplotype, haplotype_number).expect("Failed to enumerate all paths");
    let allseq = construct_sequences_from_haplotype_path(&node_info, &all_paths);
    println!("allseq: {:?}", allseq.iter().map(|(index, pathlist)| format!("index: {}, path num: {}", index, pathlist.len())).collect::<Vec<_>>().join("\n"));
    let primary_haplotypes = find_full_range_haplotypes(&node_info, &node_haplotype, &allseq, haplotype_number);
    info!("All sequences constructed: {}", primary_haplotypes.len());
    // call methylation
    let _ = call_methylation(&node_info, primary_haplotypes.clone(), output_prefix);    
    info!("Haplotype specific methylation signals exported: {}", primary_haplotypes.len());
    // write assemblies
    let output_filename = PathBuf::from(format!("{}.fasta", output_prefix.to_string_lossy()));
    let _ =write_graph_path_fasta(&primary_haplotypes, &output_filename);
    info!("All sequences written to fasta: {}", output_prefix.to_str().unwrap());
    Ok(())
}