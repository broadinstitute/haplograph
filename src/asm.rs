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
pub struct NodeInfo {
    pub seq: String,
    pub cigar: String,
    pub support_reads: String,
    pub allele_frequency: String,
}

pub fn load_graph(filename: &PathBuf) -> AnyhowResult< (HashMap<String, NodeInfo>,  HashMap<String, String>) > {
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
            };
            node_info.insert(name.to_string(), value_info);

        } else if line.starts_with('L') {
            let itemlist: Vec<&str> = line.split('\t').collect();
            let src = itemlist[1];
            let dst = itemlist[3];
            edge_info.insert(src.to_string(), dst.to_string());
        }
    }
    Ok((node_info, edge_info))
}

pub fn find_source_node(edge_info: &HashMap<String, String>) -> Vec<String> {
    let mut source_nodes = HashSet::new();
    let mut target_nodes = HashSet::new();
    for (src, dst) in edge_info.iter() {
        source_nodes.insert(src.clone());
        target_nodes.insert(dst.clone());
    }
    let start_nodes = source_nodes.difference(&target_nodes).collect::<HashSet<&String>>();
    start_nodes.into_iter().map(|s| s.clone()).collect()
}

pub fn traverse_graph(node_info: &HashMap<String, NodeInfo>, edge_info: &HashMap<String, String>) -> std::result::Result<(HashMap<Vec<String>, String>), Box<dyn Error>> {
    let source_nodes = find_source_node(edge_info);
    // DFS to traverse the graph
    let mut all_paths = Vec::new();
    for src in source_nodes {
        let mut path = Vec::new();
        path.push(src.clone());
        let mut current_node = src.clone();
        while edge_info.contains_key(&current_node) {
            let next_node = edge_info.get(&current_node).unwrap();
            path.push(next_node.clone());
            current_node = next_node.clone();
        }
        all_paths.push(path);
    }
    // get the sequence of the path
    let mut all_sequences = HashMap::new();
    for path in all_paths {
        let mut sequence = String::new();
        for node in path.iter() {
            let node_info_dict = node_info.get(node).unwrap();
            let haplotype_seq = node_info_dict.seq.clone();
            sequence += &haplotype_seq;
        }
        all_sequences.insert(path.clone(), sequence);
            
    }

    Ok(all_sequences)


}

pub fn write_graph_path_fasta(all_sequences: &HashMap<Vec<String>, String>, output_filename: &PathBuf) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;
    let chars_per_line = 60;
    for (path, sequence) in all_sequences.iter() {
        writeln!(file, ">{}", path.join("|"))?;
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

pub fn start(graph_filename: &PathBuf, output_prefix: &PathBuf) -> AnyhowResult<()> {
    let (node_info, edge_info) = load_graph(graph_filename).unwrap();
    let all_sequences = traverse_graph(&node_info, &edge_info).unwrap();
    write_graph_path_fasta(&all_sequences, output_prefix);
    Ok(())
}