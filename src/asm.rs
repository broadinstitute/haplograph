use crate::eval;
use crate::util;
use anyhow::Result as AnyhowResult;
use flate2::read::GzDecoder;
use log::{info, warn};
use ndarray::Array2;
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::path::PathBuf;

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

pub fn load_graph(
    filename: &PathBuf,
) -> AnyhowResult<(HashMap<String, NodeInfo>, HashMap<String, Vec<String>>)> {
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

            let value = json_value;
            let methyl_info = value["mod_score_dict"].as_str().unwrap_or_default();
            let methyl_info_list = if methyl_info.is_empty() {
                Vec::new()
            } else {
                methyl_info.split(",").collect::<Vec<_>>()
            };
            let methyl_info_dict = if !methyl_info_list.is_empty() {
                methyl_info_list
                    .iter()
                    .map(|x| x.split(":").collect::<Vec<_>>())
                    .collect::<Vec<_>>()
                    .iter()
                    .map(|x| (x[0].parse::<usize>().unwrap(), x[1].parse::<f32>().unwrap()))
                    .collect::<HashMap<usize, f32>>()
            } else {
                HashMap::new()
            };
            // let methyl_info_dict = methyl_info.iter().map(|x| x.split(":").collect::<Vec<_>>()).collect::<Vec<_>>().iter().map(|x| (x[0].parse::<usize>().unwrap(), x[1].parse::<f32>().unwrap())).collect::<HashMap<usize, f32>>();
            let value_info = NodeInfo {
                seq: seq.to_string(),
                cigar: value["cigar"].to_string(),
                support_reads: value["support_reads"]
                    .to_string()
                    .trim_matches('"')
                    .parse::<usize>()
                    .unwrap_or_else(|_| {
                        panic!(
                            "Support reads not found for node: {}, {}",
                            name, value["support_reads"]
                        )
                    }),
                allele_frequency: value["allele_frequency"].to_string(),
                read_names: value["read_names"].to_string(),
                methyl_info: methyl_info_dict,
            };
            node_info.insert(name.to_string(), value_info);
        } else if line.starts_with('L') {
            let itemlist: Vec<&str> = line.split('\t').collect();
            let src = itemlist[1];
            let dst = itemlist[3];
            edge_info
                .entry(src.to_string())
                .or_insert(Vec::new())
                .push(dst.to_string());
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
    let start_nodes = source_nodes
        .difference(&target_nodes)
        .collect::<HashSet<&String>>();
    start_nodes.into_iter().cloned().collect()
}

pub fn find_parallele_nodes(
    node_info: &HashMap<String, NodeInfo>,
) -> HashMap<String, HashSet<String>> {
    let mut all_nodes: HashMap<String, HashSet<String>> = HashMap::new();
    for (node, node_info) in node_info.iter() {
        let interval_name = node.split(".").collect::<Vec<_>>()[1];
        all_nodes
            .entry(interval_name.to_string().clone())
            .or_default()
            .insert(node.clone());
    }
    all_nodes
}

pub fn identify_heterozygous_nodes(
    node_info: &HashMap<String, NodeInfo>,
    hap_number: usize,
    het_fold_threshold: f64,
) -> HashMap<String, HashSet<String>> {
    let mut heterozygous_nodes: HashMap<String, HashSet<String>> = HashMap::new();
    let all_nodes = find_parallele_nodes(node_info);
    let mut interval_list = all_nodes.keys().collect::<Vec<_>>();
    interval_list.sort_by(|a, b| {
        util::split_locus(a.to_string())
            .1
            .cmp(&util::split_locus(b.to_string()).1)
    }); // ascending order a < b
    for interval_name in interval_list.iter() {
        let node_vec = all_nodes.get(interval_name.clone()).unwrap();
        if node_vec.len() < 2 {
            continue;
        }
        let node_list = node_vec.iter().map(|node| node.clone()).collect::<Vec<_>>();

        let matrix = construct_heterozygous_nodes_matrix(node_info, node_list.clone());
        let node_list_filtered = util::permutation_test(&matrix, 0.1, 50, node_list.clone());

        if node_list_filtered.len() == hap_number {
            heterozygous_nodes
                .entry(interval_name.to_string().clone())
                .or_default()
                .extend(node_list_filtered.iter().map(|node| node.clone()));
        }
    }
    heterozygous_nodes
}

fn hungarian_minimum_assignment(cost: &[Vec<i64>]) -> Vec<usize> {
    let size = cost.len();
    if size == 0 {
        return Vec::new();
    }

    let inf = i64::MAX / 4;
    let mut u = vec![0_i64; size + 1];
    let mut v = vec![0_i64; size + 1];
    let mut p = vec![0_usize; size + 1];
    let mut way = vec![0_usize; size + 1];

    for i in 1..=size {
        p[0] = i;
        let mut j0 = 0_usize;
        let mut minv = vec![inf; size + 1];
        let mut used = vec![false; size + 1];

        loop {
            used[j0] = true;
            let i0 = p[j0];
            let mut delta = inf;
            let mut j1 = 0_usize;

            for j in 1..=size {
                if used[j] {
                    continue;
                }
                let cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                if cur < minv[j] {
                    minv[j] = cur;
                    way[j] = j0;
                }
                if minv[j] < delta {
                    delta = minv[j];
                    j1 = j;
                }
            }

            for j in 0..=size {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }

            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }

        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }

    let mut assignment = vec![0_usize; size];
    for j in 1..=size {
        if p[j] != 0 {
            assignment[p[j] - 1] = j - 1;
        }
    }
    assignment
}

fn hungarian_maximum_matching(weights: &[Vec<i64>]) -> Vec<Option<usize>> {
    if weights.is_empty() {
        return Vec::new();
    }

    let row_count = weights.len();
    let col_count = weights.iter().map(|row| row.len()).max().unwrap_or(0);
    if col_count == 0 {
        return vec![None; row_count];
    }

    let size = row_count.max(col_count);
    let max_weight = weights
        .iter()
        .flat_map(|row| row.iter())
        .copied()
        .max()
        .unwrap_or(0);
    let mut cost = vec![vec![max_weight; size]; size];

    for (row_index, row) in weights.iter().enumerate() {
        for (col_index, weight) in row.iter().enumerate() {
            cost[row_index][col_index] = max_weight - *weight;
        }
    }

    hungarian_minimum_assignment(&cost)
        .into_iter()
        .take(row_count)
        .map(|col_index| {
            if col_index < col_count {
                Some(col_index)
            } else {
                None
            }
        })
        .collect()
}

fn maximum_overlap_matching(
    haplotype_read_sets: &[(usize, HashSet<String>)],
    candidate_node_reads: &[(String, HashSet<String>)],
) -> HashMap<usize, String> {
    if haplotype_read_sets.is_empty() || candidate_node_reads.is_empty() {
        return HashMap::new();
    }

    let weights = haplotype_read_sets
        .iter()
        .map(|(_, hap_reads)| {
            candidate_node_reads
                .iter()
                .map(|(_, node_reads)| hap_reads.intersection(node_reads).count() as i64)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let assignments = hungarian_maximum_matching(&weights);
    let mut matched_nodes = HashMap::new();
    for (row_index, assignment) in assignments.into_iter().enumerate() {
        if let Some(col_index) = assignment {
            if weights[row_index]
                .get(col_index)
                .copied()
                .unwrap_or_default()
                > 0
            {
                matched_nodes.insert(
                    haplotype_read_sets[row_index].0,
                    candidate_node_reads[col_index].0.clone(),
                );
            }
        }
    }

    matched_nodes
}

fn supported_haplotypes_for_node(
    haplotype_read_sets: &[(usize, HashSet<String>)],
    node_reads: &HashSet<String>,
) -> Vec<usize> {
    let mut supported_haplotypes = haplotype_read_sets
        .iter()
        .filter_map(|(haplotype, haplotype_reads)| {
            haplotype_reads
                .intersection(node_reads)
                .next()
                .map(|_| *haplotype)
        })
        .collect::<Vec<_>>();
    supported_haplotypes.sort_unstable();
    supported_haplotypes
}

fn interval_from_node_id(node_id: &str) -> String {
    node_id
        .split('.')
        .nth(1)
        .unwrap_or(node_id)
        .to_string()
}

fn haplotype_read_overlap(
    node_id: &str,
    haplotype: usize,
    node_info: &HashMap<String, NodeInfo>,
    haplotype_reads: &HashMap<usize, HashSet<String>>,
) -> usize {
    let node_reads = get_read_name_list(node_info, node_id.to_string());
    haplotype_reads
        .get(&haplotype)
        .map(|hap_reads| node_reads.intersection(hap_reads).count())
        .unwrap_or(0)
}

/// Per interval, each haplotype index may label at most one node.
/// A node may still carry multiple haplotype labels (e.g. shared/deletion windows).
/// When the same haplotype labels multiple nodes, keep the node most similar to that
/// haplotype's known read set and remove the label from the others.
pub fn deduplicate_interval_node_haplotype_labels(
    node_haplotype: HashMap<String, HashSet<usize>>,
    node_info: &HashMap<String, NodeInfo>,
    haplotype_reads: &HashMap<usize, HashSet<String>>,
    hap_number: usize,
) -> HashMap<String, HashSet<usize>> {
    let mut intervals: HashMap<String, Vec<String>> = HashMap::new();
    for node_id in node_haplotype.keys() {
        intervals
            .entry(interval_from_node_id(node_id))
            .or_default()
            .push(node_id.clone());
    }

    let mut deduped = node_haplotype;

    for (interval, nodes) in intervals {
        // Same haplotype label on multiple nodes in this interval -> keep the best-supported node.
        for hap in 0..hap_number {
            let mut claiming_nodes: Vec<String> = nodes
                .iter()
                .filter(|node_id| {
                    deduped
                        .get(*node_id)
                        .is_some_and(|labels| labels.contains(&hap))
                })
                .cloned()
                .collect();
            if claiming_nodes.len() <= 1 {
                continue;
            }

            claiming_nodes.sort_by(|a, b| {
                haplotype_read_overlap(b, hap, node_info, haplotype_reads)
                    .cmp(&haplotype_read_overlap(a, hap, node_info, haplotype_reads))
                    .then(a.cmp(b))
            });
            let keep_node = claiming_nodes[0].clone();
            for node_id in claiming_nodes.iter().skip(1) {
                if let Some(labels) = deduped.get_mut(node_id) {
                    labels.remove(&hap);
                    if labels.is_empty() {
                        deduped.remove(node_id);
                    }
                }
            }
            info!(
                "interval {interval}: haplotype {hap} kept on {keep_node}, removed from {} other node(s)",
                claiming_nodes.len() - 1
            );
        }
    }

    deduped
}

pub fn log_interval_duplicate_haplotype_labels(
    node_haplotype: &HashMap<String, HashSet<usize>>,
    hap_number: usize,
    stage: &str,
) {
    let mut intervals: HashMap<String, Vec<(String, HashSet<usize>)>> = HashMap::new();
    for (node_id, labels) in node_haplotype {
        intervals
            .entry(interval_from_node_id(node_id))
            .or_default()
            .push((node_id.clone(), labels.clone()));
    }

    for (interval, nodes) in intervals {
        let mut hap_to_nodes: HashMap<usize, Vec<String>> = HashMap::new();
        for (node_id, labels) in &nodes {
            for hap in labels {
                hap_to_nodes.entry(*hap).or_default().push(node_id.clone());
            }
        }

        for hap in 0..hap_number {
            if let Some(node_ids) = hap_to_nodes.get(&hap) {
                if node_ids.len() > 1 {
                    warn!(
                        "{stage} interval {interval}: haplotype {hap} on {} nodes: {}",
                        node_ids.len(),
                        node_ids.join(", ")
                    );
                } else if node_ids.is_empty() {
                    warn!("{stage} interval {interval}: haplotype {hap} has no node");
                }
            } else {
                warn!("{stage} interval {interval}: haplotype {hap} has no node");
            }
        }
    }
}

pub fn assign_haplotype_reads(
    node_info: &HashMap<String, NodeInfo>,
    heterozygous_nodes: &HashMap<String, HashSet<String>>,
    hap_number: usize,
) -> HashMap<usize, HashSet<String>> {
    // assign reads to haplotypes
    let mut haplotype_reads = HashMap::new();
    let mut interval_list = heterozygous_nodes.keys().collect::<Vec<_>>();
    interval_list.sort_by(|a, b| {
        util::split_locus(a.to_string())
            .1
            .cmp(&util::split_locus(b.to_string()).1)
    });
    for interval_name in interval_list.iter() {
        let mut node_set = heterozygous_nodes
            .get(interval_name.clone())
            .unwrap()
            .clone()
            .iter()
            .cloned()
            .collect::<Vec<_>>();
        let mut node_read_sets = node_set
            .drain(..)
            .map(|node| {
                let read_name_list = get_read_name_list(node_info, node.clone());
                (node, read_name_list)
            })
            .collect::<Vec<_>>();
        node_read_sets.sort_by(|a, b| a.1.len().cmp(&b.1.len()).then(a.0.cmp(&b.0)));

        if haplotype_reads.is_empty() {
            if node_read_sets.len() < hap_number {
                continue;
            }
            for (index, (_, read_name_list)) in node_read_sets.iter().take(hap_number).enumerate() {
                haplotype_reads
                    .entry(index)
                    .or_insert(HashSet::new())
                    .extend(read_name_list.clone());
            }
        } else {
            let mut haplotype_order = haplotype_reads.keys().copied().collect::<Vec<_>>();
            haplotype_order.sort_unstable();
            let haplotype_read_sets = haplotype_order
                .iter()
                .filter_map(|hap| haplotype_reads.get(hap).cloned().map(|reads| (*hap, reads)))
                .collect::<Vec<_>>();
            let matched_nodes = maximum_overlap_matching(&haplotype_read_sets, &node_read_sets);

            for (haplotype, node_id) in matched_nodes {
                if let Some((_, read_name_list)) =
                    node_read_sets.iter().find(|(node, _)| *node == node_id)
                {
                    haplotype_reads
                        .entry(haplotype)
                        .or_insert(HashSet::new())
                        .extend(read_name_list.clone());
                }
            }
        }
    }

    haplotype_reads
}

pub fn assign_node_to_reads(
    node_info: &HashMap<String, NodeInfo>,
) -> HashMap<String, HashSet<String>> {
    let mut read_to_nodes = HashMap::new();
    let nodelist = node_info.keys().collect::<Vec<_>>();
    for node in nodelist {
        let read_names = get_read_name_list(node_info, node.clone());
        for read in read_names {
            read_to_nodes
                .entry(read)
                .or_insert(HashSet::new())
                .insert(node.clone());
        }
    }
    read_to_nodes
}

pub fn assign_haplotype_nodes(
    node_info: &HashMap<String, NodeInfo>,
    haplotype_reads: &HashMap<usize, HashSet<String>>,
) -> HashMap<usize, HashSet<String>> {
    let read_to_nodes = assign_node_to_reads(node_info);
    let mut haplotype_nodes = HashMap::new();
    for (haplotype, reads) in haplotype_reads.iter() {
        for r in reads {
            haplotype_nodes
                .entry(*haplotype)
                .or_insert(HashSet::new())
                .extend(read_to_nodes.get(r).unwrap().iter().cloned());
        }
    }
    haplotype_nodes
}

pub fn find_unassigned_reads(
    node_info: &HashMap<String, NodeInfo>,
    haplotype_reads: &HashMap<usize, HashSet<String>>,
) -> HashSet<String> {
    let mut assigned_reads = HashSet::new();
    for reads in haplotype_reads.values() {
        assigned_reads.extend(reads.iter().cloned());
    }
    let all_reads = find_all_reads(node_info);
    all_reads.difference(&assigned_reads).cloned().collect()
}

pub fn construct_heterozygous_nodes_matrix(
    node_info: &HashMap<String, NodeInfo>,
    node_list: Vec<String>,
) -> Array2<f64> {
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

pub fn filter_heterozygous_nodes(
    node_info: &HashMap<String, NodeInfo>,
    heterozygous_nodes: &HashMap<String, HashSet<String>>,
) -> HashMap<String, HashSet<String>> {
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
        filtered_heterozygous_nodes
            .entry(interval_name.clone())
            .or_insert(HashSet::new())
            .insert(node_id.clone());
    }
    filtered_heterozygous_nodes
}

pub fn assign_unassigned_reads(
    node_info: &HashMap<String, NodeInfo>,
    haplotype_reads: &HashMap<usize, HashSet<String>>,
) -> (
    HashMap<usize, HashSet<String>>,
    HashMap<usize, HashSet<String>>,
) {
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
            let overlap_nodes = hap_nodes
                .intersection(read_nodes)
                .cloned()
                .collect::<HashSet<_>>();
            if !overlap_nodes.is_empty() {
                haplotype_reads_new
                    .entry(*haplotype)
                    .or_default()
                    .insert(read.clone());
                haplotype_nodes_new
                    .entry(*haplotype)
                    .or_default()
                    .extend(read_nodes.clone());
                find_haplotype = true;
            }
        }
        // homologous reads
        if !find_haplotype {
            for hap in haplotype_reads_new.clone().keys() {
                haplotype_reads_new
                    .entry(*hap)
                    .or_default()
                    .insert(read.clone());
                haplotype_nodes_new
                    .entry(*hap)
                    .or_default()
                    .extend(read_nodes.clone());
            }
        }
    }
    (haplotype_reads_new, haplotype_nodes_new)
}

pub fn find_most_supported_path(
    node_info: &HashMap<String, NodeInfo>,
) -> HashMap<String, HashSet<usize>> {
    let mut node_haplotype = HashMap::new();
    let all_nodes = find_parallele_nodes(node_info);
    let mut interval_list = all_nodes.keys().collect::<Vec<_>>();
    interval_list.sort(); // ascending order a < b
                          // should also consider the edge connectivity
    let haplotype = 0;
    for (interval_name, node_vec) in all_nodes.iter() {
        if node_vec.len() < 2 {
            let node_id = node_vec.iter().next().unwrap().clone();
            node_haplotype
                .entry(node_id)
                .or_insert(HashSet::new())
                .insert(haplotype);
        } else {
            let mut node_support = Vec::new();
            for node in node_vec.iter() {
                let read_name_list = get_read_name_list(node_info, node.clone());
                node_support.push((node.clone(), read_name_list.len()));
            }
            node_support.sort_by(|a, b| b.1.cmp(&a.1)); // descending order a > b
            let node_id = node_support[0].0.clone();
            node_haplotype
                .entry(node_id)
                .or_insert(HashSet::new())
                .insert(haplotype);
        }
    }

    node_haplotype
}

pub fn assign_haplotype_to_nodes(
    haplotype_nodes: &HashMap<usize, HashSet<String>>,
) -> HashMap<String, HashSet<usize>> {
    let mut node_haplotype = HashMap::new();
    for (haplotype, nodes) in haplotype_nodes.iter() {
        for node in nodes {
            node_haplotype
                .entry(node.clone())
                .or_insert(HashSet::new())
                .insert(*haplotype);
        }
    }
    node_haplotype
}

pub fn find_all_reads(node_info: &HashMap<String, NodeInfo>) -> HashSet<String> {
    let mut all_reads = HashSet::new();
    for (node, node_infomation) in node_info.iter() {
        let read_names = get_read_name_list(node_info, node.clone());
        all_reads.extend(read_names.iter().cloned());
    }
    all_reads
}

pub fn get_read_name_list(
    node_info: &HashMap<String, NodeInfo>,
    node_id: String,
) -> HashSet<String> {
    let read_names = node_info.get(&node_id).unwrap().read_names.clone();
    let read_names_list = read_names
        .split(",")
        .collect::<Vec<_>>()
        .iter()
        .map(|x| x.to_string())
        .collect::<HashSet<_>>();
    let read_names_list_clone = read_names_list
        .iter()
        .map(|x| {
            x.split("|").collect::<Vec<_>>()[0]
                .to_string()
                .replace("\"", "")
                .clone()
        })
        .collect::<HashSet<_>>();
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
        let mut haplotype_intersection_src =
            constrain_haplotype_index(&haplotype_index, node_haplotype.get(&src));
        if !haplotype_intersection_src.is_empty() {
            let mut path = Vec::new();
            path.push(src.clone());
            dfs_traverse_with_haplotype_constrains(
                &src,
                edge_info,
                &mut path,
                &mut all_paths,
                &mut haplotype_intersection_src,
                node_haplotype,
            );
        }
    }

    Ok(all_paths)
}

fn constrain_haplotype_index(
    haplotype_index: &HashSet<usize>,
    node_haplotype: Option<&HashSet<usize>>,
) -> HashSet<usize> {
    match node_haplotype {
        Some(node_haplotypes) => haplotype_index
            .intersection(node_haplotypes)
            .cloned()
            .collect::<HashSet<_>>(),
        None => haplotype_index.clone(),
    }
}

// /// Recursive DFS to find all paths from a starting node
fn dfs_traverse_with_haplotype_constrains(
    current_node: &String,
    edge_info: &HashMap<String, Vec<String>>,
    current_path: &mut Vec<String>,
    all_paths: &mut Vec<(Vec<String>, HashSet<usize>)>,
    haplotype_index: &mut HashSet<usize>,
    node_haplotype: &HashMap<String, HashSet<usize>>,
) {
    if !edge_info.contains_key(current_node) {
        if !haplotype_index.is_empty() {
            all_paths.push((current_path.clone(), haplotype_index.clone()));
        }
        return;
    }

    if haplotype_index.is_empty() {
        return;
    }

    let next_nodes = edge_info.get(current_node).unwrap();
    for next_node in next_nodes {
        let mut haplotype_intersection_clone =
            constrain_haplotype_index(haplotype_index, node_haplotype.get(next_node));
        current_path.push(next_node.clone());
        dfs_traverse_with_haplotype_constrains(
            next_node,
            edge_info,
            current_path,
            all_paths,
            &mut haplotype_intersection_clone,
            node_haplotype,
        );
        current_path.pop(); // Backtrack
    }
}

pub fn construct_sequences_from_haplotype_path(
    node_info: &HashMap<String, NodeInfo>,
    all_paths: &Vec<(Vec<String>, HashSet<usize>)>,
) -> HashMap<usize, Vec<(Vec<String>, String, HashSet<String>)>> {
    // Get the sequence of each path
    let mut all_sequences = HashMap::new();
    for (path_index, (path, haplotype_index)) in all_paths.iter().enumerate() {
        let mut sequence = String::new();
        let mut read_names = HashSet::new();

        for node in path.iter() {
            let node_info_dict = node_info.get(node).unwrap_or_else(|| {
                panic!(
                    "Node {} not found in node_info, path Number: {:?}",
                    node, path
                )
            });
            let read_names_list_clone = get_read_name_list(node_info, node.clone());
            read_names.extend(read_names_list_clone.clone());
            let haplotype_seq = node_info_dict.seq.clone();
            sequence += &haplotype_seq;
            // total_supported_reads += supported_reads;
        }
        for hap_ind in haplotype_index.iter() {
            info!(
                "path: {:?}, haplotype_index: {:?}, sequence length: {:?}",
                path.len(),
                hap_ind,
                sequence.len()
            );
            all_sequences.entry(*hap_ind).or_insert(Vec::new()).push((
                path.clone(),
                sequence.clone(),
                read_names.clone(),
            ));
        }
    }
    all_sequences
}

pub fn write_graph_path_fasta(
    all_sequences: &HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)>,
    output_filename: &PathBuf,
) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;
    let chars_per_line = 60;
    for (index, (path, sequence, supported_reads, supports, span)) in all_sequences.iter() {
        let chromosome = path[0].split(".").collect::<Vec<_>>()[1]
            .split(":")
            .collect::<Vec<_>>()[0];
        let (start, end) =
            eval::find_alignment_intervals(path.iter().map(|x| x.as_str()).collect::<Vec<_>>())?;
        writeln!(
            file,
            ">{}:{}-{}.{}\tSupports:{}\t{}",
            chromosome,
            start,
            end,
            index,
            supported_reads.len(),
            path.join("|")
        )?;

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

pub fn get_supports(node_info: &HashMap<String, NodeInfo>, path: &Vec<String>) -> usize {
    let mut supports = 0_usize;
    for node in path.iter() {
        let read_set = get_read_name_list(node_info, node.clone());
        supports += read_set.len();
    }
    supports
}

pub fn find_full_range_haplotypes(
    node_info: &HashMap<String, NodeInfo>,
    node_haplotype: &HashMap<String, HashSet<usize>>,
    all_sequences: &HashMap<usize, Vec<(Vec<String>, String, HashSet<String>)>>,
) -> HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)> {
    // sort all_sequences by the spanning length and the supported_reads
    let node_list = node_info.keys().collect::<Vec<_>>();
    // let full_range = eval::find_alignment_intervals(node_list.iter().map(|x| x.as_str()).collect::<Vec<_>>()).unwrap();
    let mut best_paths = HashMap::new();
    for (hap_index, path_list) in all_sequences.iter() {
        // for each haplotype, select the best path
        let mut full_sequences = Vec::new();
        for (index, (path, sequence, supported_reads)) in path_list.iter().enumerate() {
            let (start, end) =
                eval::find_alignment_intervals(path.iter().map(|x| x.as_str()).collect::<Vec<_>>())
                    .unwrap();
            let supports = get_supports(node_info, path);
            info!(
                "path: {:?}, supports: {:?}, supported_reads: {:?}, span: {:?}",
                hap_index,
                supports,
                supported_reads.len(),
                end - start
            );

            full_sequences.push((
                path.clone(),
                sequence.clone(),
                supported_reads.clone(),
                supports,
                end - start,
            ));
        }
        //first compare the 4th element, then the 3th element
        // full_sequences.sort_by(|a, b| b.3.cmp(&a.3).then(b.4.cmp(&a.4)));
        full_sequences.sort_by(|a, b| b.2.len().cmp(&a.2.len()).then(b.3.cmp(&a.3)));
        best_paths.insert(*hap_index, full_sequences[0].clone());
    }
    best_paths
}

pub fn find_parallele_nodes_from_nodelist(
    nodelist: &Vec<String>,
) -> HashMap<String, HashSet<String>> {
    let mut interval_node = HashMap::new();
    for node in nodelist {
        let interval = node.split(".").collect::<Vec<_>>()[1];
        interval_node
            .entry(interval.to_string())
            .or_insert(HashSet::new())
            .insert(node.clone());
    }
    interval_node
}

pub fn filter_haplotype_nodes(
    node_info: &HashMap<String, NodeInfo>,
    haplotype_nodes: &HashMap<usize, HashSet<String>>,
    haplotype_reads: &HashMap<usize, HashSet<String>>,
) -> HashMap<usize, HashSet<String>> {
    let mut filtered_haplotype_nodes = HashMap::new();

    let all_nodes = haplotype_nodes
        .values()
        .flat_map(|nodes| nodes.iter().cloned())
        .collect::<Vec<_>>();
    let interval_node = find_parallele_nodes_from_nodelist(&all_nodes);
    let mut interval_list = interval_node.keys().collect::<Vec<_>>();
    interval_list.sort_by(|a, b| {
        util::split_locus(a.to_string())
            .1
            .cmp(&util::split_locus(b.to_string()).1)
    });

    let mut haplotype_order = haplotype_reads.keys().copied().collect::<Vec<_>>();
    haplotype_order.sort_unstable();
    let haplotype_read_sets = haplotype_order
        .iter()
        .filter_map(|hap| haplotype_reads.get(hap).cloned().map(|reads| (*hap, reads)))
        .collect::<Vec<_>>();

    for interval in interval_list.iter() {
        let nodes = interval_node.get(*interval).unwrap().clone();
        let mut node_read_sets = nodes
            .iter()
            .cloned()
            .map(|node| {
                let node_read_names = get_read_name_list(node_info, node.clone());
                (node, node_read_names)
            })
            .collect::<Vec<_>>();
        node_read_sets.sort_by(|a, b| a.1.len().cmp(&b.1.len()).then(a.0.cmp(&b.0)));

        if node_read_sets.len() == 1 {
            let (node_id, node_reads) = &node_read_sets[0];
            let supported_haplotypes =
                supported_haplotypes_for_node(&haplotype_read_sets, node_reads);
            if !supported_haplotypes.is_empty() {
                // Shared deletion-only windows often collapse to a single node.
                // Keep that node on every supported haplotype rather than forcing
                // a one-to-one assignment that would drop one haplotype path.
                for haplotype in supported_haplotypes {
                    filtered_haplotype_nodes
                        .entry(haplotype)
                        .or_insert(HashSet::new())
                        .insert(node_id.clone());
                }
                continue;
            }
        }

        let matched_nodes = maximum_overlap_matching(&haplotype_read_sets, &node_read_sets);
        if matched_nodes.is_empty() {
            warn!(
                "interval: {}, nodes: {:?}",
                interval,
                nodes.iter().cloned().collect::<Vec<_>>().join(", ")
            );
            continue;
        }

        for (haplotype, node_id) in matched_nodes {
            filtered_haplotype_nodes
                .entry(haplotype)
                .or_insert(HashSet::new())
                .insert(node_id);
        }
    }
    filtered_haplotype_nodes
}

pub fn find_node_haplotype(
    node_info: &HashMap<String, NodeInfo>,
    hap_number: usize,
    het_fold_threshold: f64,
) -> (
    HashMap<usize, HashSet<String>>,
    HashMap<String, HashSet<usize>>,
) {
    if hap_number == 1 {
        let node_haplotype = find_most_supported_path(node_info);
        return (HashMap::new(), node_haplotype);
    }
    if hap_number < 1 {
        warn!("Unsupported hap_number: {}", hap_number);
        return (HashMap::new(), HashMap::new());
    }

    let heterozygous_nodes = identify_heterozygous_nodes(node_info, hap_number, het_fold_threshold);
    info!("heterozygous_nodes: {:?}", heterozygous_nodes.len());
    let filtered_heterozygous_nodes = filter_heterozygous_nodes(node_info, &heterozygous_nodes);
    info!(
        "filtered_heterozygous_nodes: {:?}",
        filtered_heterozygous_nodes.len()
    );

    let haplotype_reads =
        assign_haplotype_reads(node_info, &filtered_heterozygous_nodes, hap_number);

    let (haplotype_reads_new, haplotype_nodes_new) =
        assign_unassigned_reads(node_info, &haplotype_reads);
    info!(
        "haplotype_reads: {:?}",
        haplotype_reads
            .iter()
            .map(|(hap, reads)| format!("hap: {}, reads: {}", hap, reads.len()))
            .collect::<Vec<_>>()
            .join(", ")
    );
    info!(
        "haplotype_reads_new: {:?}",
        haplotype_reads_new
            .iter()
            .map(|(hap, reads)| format!("hap: {}, reads: {}", hap, reads.len()))
            .collect::<Vec<_>>()
            .join(", ")
    );

    if haplotype_reads_new.is_empty() {
        let node_haplotype = find_most_supported_path(node_info);
        (haplotype_reads_new, node_haplotype)
    } else {
        let filtered_haplotype_nodes =
            filter_haplotype_nodes(node_info, &haplotype_nodes_new, &haplotype_reads_new);
        let mut node_haplotype = assign_haplotype_to_nodes(&filtered_haplotype_nodes);
        log_interval_duplicate_haplotype_labels(&node_haplotype, hap_number, "before dedup");
        node_haplotype = deduplicate_interval_node_haplotype_labels(
            node_haplotype,
            node_info,
            &haplotype_reads_new,
            hap_number,
        );
        log_interval_duplicate_haplotype_labels(&node_haplotype, hap_number, "after dedup");
        (haplotype_reads_new, node_haplotype)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_node(reads: &[&str]) -> NodeInfo {
        NodeInfo {
            seq: "A".to_string(),
            cigar: "1=".to_string(),
            support_reads: reads.len(),
            allele_frequency: String::new(),
            read_names: reads
                .iter()
                .map(|read| format!("{}|0", read))
                .collect::<Vec<_>>()
                .join(","),
            methyl_info: HashMap::new(),
        }
    }

    #[test]
    fn hungarian_matching_finds_global_optimum() {
        let weights = vec![vec![10, 9], vec![8, 1]];
        let assignment = hungarian_maximum_matching(&weights);
        assert_eq!(assignment, vec![Some(1), Some(0)]);
    }

    #[test]
    fn assign_haplotype_reads_uses_unique_interval_matching() {
        let mut node_info = HashMap::new();
        node_info.insert(
            "graph.chr1:0-10.a".to_string(),
            test_node(&["u", "v", "w", "x"]),
        );
        node_info.insert(
            "graph.chr1:0-10.b".to_string(),
            test_node(&["a", "b", "c", "d", "e"]),
        );
        node_info.insert(
            "graph.chr1:10-20.a".to_string(),
            test_node(&["a", "b", "c", "d"]),
        );
        node_info.insert(
            "graph.chr1:10-20.b".to_string(),
            test_node(&["a", "b", "c", "d", "e", "u", "v", "w", "x"]),
        );

        let heterozygous_nodes = HashMap::from([
            (
                "chr1:0-10".to_string(),
                HashSet::from([
                    "graph.chr1:0-10.a".to_string(),
                    "graph.chr1:0-10.b".to_string(),
                ]),
            ),
            (
                "chr1:10-20".to_string(),
                HashSet::from([
                    "graph.chr1:10-20.a".to_string(),
                    "graph.chr1:10-20.b".to_string(),
                ]),
            ),
        ]);

        let haplotype_reads = assign_haplotype_reads(&node_info, &heterozygous_nodes, 2);
        assert!(haplotype_reads.get(&0).unwrap().contains("a"));
    }

    #[test]
    fn filter_haplotype_nodes_uses_unique_interval_matching() {
        let mut node_info = HashMap::new();
        node_info.insert(
            "graph.chr1:10-20.a".to_string(),
            test_node(&["a", "b", "c", "d"]),
        );
        node_info.insert(
            "graph.chr1:10-20.b".to_string(),
            test_node(&["a", "b", "c", "d", "e", "u", "v", "w", "x"]),
        );

        let haplotype_nodes = HashMap::from([
            (
                0,
                HashSet::from([
                    "graph.chr1:10-20.a".to_string(),
                    "graph.chr1:10-20.b".to_string(),
                ]),
            ),
            (
                1,
                HashSet::from([
                    "graph.chr1:10-20.a".to_string(),
                    "graph.chr1:10-20.b".to_string(),
                ]),
            ),
        ]);
        let haplotype_reads = HashMap::from([
            (
                0,
                HashSet::from(["a", "b", "c", "d", "e"].map(String::from)),
            ),
            (1, HashSet::from(["u", "v", "w", "x"].map(String::from))),
        ]);

        let filtered = filter_haplotype_nodes(&node_info, &haplotype_nodes, &haplotype_reads);
        assert_eq!(
            filtered.get(&0).unwrap(),
            &HashSet::from(["graph.chr1:10-20.a".to_string()])
        );
        assert_eq!(
            filtered.get(&1).unwrap(),
            &HashSet::from(["graph.chr1:10-20.b".to_string()])
        );
    }

    #[test]
    fn deduplicate_interval_node_haplotype_labels_allows_multi_label_node() {
        let mut node_info = HashMap::new();
        node_info.insert(
            "graph.chr1:10-20.a".to_string(),
            test_node(&["a", "b", "c", "d"]),
        );
        node_info.insert(
            "graph.chr1:10-20.b".to_string(),
            test_node(&["a", "b", "c", "d", "e", "u", "v", "w", "x"]),
        );
        let haplotype_reads = HashMap::from([
            (
                0,
                HashSet::from(["a", "b", "c", "d", "e"].map(String::from)),
            ),
            (1, HashSet::from(["u", "v", "w", "x"].map(String::from))),
        ]);
        let node_haplotype = HashMap::from([
            ("graph.chr1:10-20.a".to_string(), HashSet::from([0])),
            ("graph.chr1:10-20.b".to_string(), HashSet::from([0, 1])),
        ]);

        let deduped = deduplicate_interval_node_haplotype_labels(
            node_haplotype,
            &node_info,
            &haplotype_reads,
            2,
        );
        assert_eq!(
            deduped.get("graph.chr1:10-20.a").unwrap(),
            &HashSet::from([0])
        );
        assert_eq!(
            deduped.get("graph.chr1:10-20.b").unwrap(),
            &HashSet::from([0, 1])
        );
    }

    #[test]
    fn deduplicate_interval_node_haplotype_labels_keeps_best_node_per_haplotype() {
        let mut node_info = HashMap::new();
        node_info.insert(
            "graph.chr1:0-10.a".to_string(),
            test_node(&["u", "v"]),
        );
        node_info.insert(
            "graph.chr1:0-10.b".to_string(),
            test_node(&["a", "b", "c", "d", "e"]),
        );
        let haplotype_reads = HashMap::from([(
            0,
            HashSet::from(["a", "b", "c", "d", "e"].map(String::from)),
        )]);
        let node_haplotype = HashMap::from([
            ("graph.chr1:0-10.a".to_string(), HashSet::from([0])),
            ("graph.chr1:0-10.b".to_string(), HashSet::from([0])),
        ]);

        let deduped = deduplicate_interval_node_haplotype_labels(
            node_haplotype,
            &node_info,
            &haplotype_reads,
            2,
        );
        assert_eq!(deduped.len(), 1);
        assert_eq!(
            deduped.get("graph.chr1:0-10.b").unwrap(),
            &HashSet::from([0])
        );
    }

    #[test]
    fn filter_haplotype_nodes_keeps_shared_singleton_node_on_multiple_haplotypes() {
        let mut node_info = HashMap::new();
        node_info.insert(
            "graph.chr1:20-30.shared".to_string(),
            test_node(&["a", "b", "u", "v"]),
        );

        let haplotype_nodes = HashMap::from([
            (0, HashSet::from(["graph.chr1:20-30.shared".to_string()])),
            (1, HashSet::from(["graph.chr1:20-30.shared".to_string()])),
        ]);
        let haplotype_reads = HashMap::from([
            (0, HashSet::from(["a", "b", "c"].map(String::from))),
            (1, HashSet::from(["u", "v", "w"].map(String::from))),
        ]);

        let filtered = filter_haplotype_nodes(&node_info, &haplotype_nodes, &haplotype_reads);
        let expected = HashSet::from(["graph.chr1:20-30.shared".to_string()]);
        assert_eq!(filtered.get(&0).unwrap(), &expected);
        assert_eq!(filtered.get(&1).unwrap(), &expected);
    }

    #[test]
    fn path_enumeration_keeps_unlabeled_intermediate_nodes() {
        let edge_info = HashMap::from([
            ("n0".to_string(), vec!["n1".to_string()]),
            ("n1".to_string(), vec!["n2".to_string()]),
        ]);
        let node_haplotype = HashMap::from([
            ("n0".to_string(), HashSet::from([0])),
            ("n2".to_string(), HashSet::from([0])),
        ]);

        let all_paths =
            enumerate_all_paths_with_haplotype(&HashMap::new(), &edge_info, &node_haplotype, 1)
                .unwrap();

        assert_eq!(
            all_paths,
            vec![(
                vec!["n0".to_string(), "n1".to_string(), "n2".to_string()],
                HashSet::from([0])
            )]
        );
    }

    #[test]
    fn path_enumeration_keeps_unlabeled_source_nodes() {
        let edge_info = HashMap::from([("n0".to_string(), vec!["n1".to_string()])]);
        let node_haplotype = HashMap::from([("n1".to_string(), HashSet::from([0]))]);

        let all_paths =
            enumerate_all_paths_with_haplotype(&HashMap::new(), &edge_info, &node_haplotype, 1)
                .unwrap();

        assert_eq!(
            all_paths,
            vec![(vec!["n0".to_string(), "n1".to_string()], HashSet::from([0]))]
        );
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
        if c.is_ascii_digit() {
            num.push(c);
        } else if !num.is_empty() {
            let length = num.parse::<usize>().unwrap();
            operations.push((length, c));
            num.clear();
        }
    }

    // Process each operation
    let mut position_mapping = HashMap::new();

    for (length, op) in operations {
        match op {
            '=' | 'M' => {
                // Match
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    position_mapping.insert(alt_pos + i, pos);
                }
                ref_pos += length;
                alt_pos += length;
            }
            'X' => {
                // Mismatch
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    position_mapping.insert(alt_pos + i, pos);
                }
                ref_pos += length;
                alt_pos += length;
            }
            'I' => {
                // Insertion
                let pos = ref_start + ref_pos - 1; // Mapping to one base pair before
                for i in 0..length {
                    position_mapping.insert(alt_pos + i, pos);
                }
                alt_pos += length;
            }
            'D' => {
                // Deletion, nothing will map back to reference coordinates
                ref_pos += length;
            }
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
        methyl_info_vec.push((
            chromosome,
            ref_pos_start,
            ref_pos_end,
            asm_pos_start,
            asm_pos_end,
            motif,
            score,
            coverage,
        ));
    }
    methyl_info_vec.sort_by_key(|item| item.1);
    for methyl_list in methyl_info_vec.iter() {
        let (
            chromosome,
            ref_pos_start,
            ref_pos_end,
            asm_pos_start,
            asm_pos_end,
            motif,
            methyl_rate,
            coverage,
        ) = methyl_list;
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chromosome,
            ref_pos_start,
            ref_pos_end,
            methyl_rate,
            asm_pos_start,
            asm_pos_end,
            motif,
            coverage
        )?;
    }

    Ok(())
}

pub fn call_methylation(
    node_info: &HashMap<String, NodeInfo>,
    all_paths: HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)>,
    output_prefix: &PathBuf,
) {
    // call methylation
    for (hap_index, (path, sequence, supported_reads, supports, span)) in all_paths.iter() {
        let mut methyl_info_dict = HashMap::new();
        let mut spos = 0;
        let chromosome = path[0].split(".").collect::<Vec<_>>()[1]
            .split(":")
            .collect::<Vec<_>>()[0];
        for node in path.iter() {
            let methyl_info = node_info.get(node).unwrap().methyl_info.clone();
            let cigar = node_info.get(node).unwrap().cigar.clone();
            let ref_start = eval::find_alignment_intervals(
                [node.clone()]
                    .iter()
                    .map(|x| x.as_str())
                    .collect::<Vec<_>>(),
            )
            .unwrap()
            .0;
            let position_mapping = mapping_to_reference_coordinates(&cigar, ref_start);
            // println!("position_mapping: {:?}", position_mapping);
            let node_seq = node_info.get(node).unwrap().seq.clone();
            let node_coverage = node_info.get(node).unwrap().support_reads;
            for (pos, score) in methyl_info.into_iter() {
                let asm_pos = pos + spos;
                let ref_pos = *position_mapping.get(&pos).unwrap_or(&0);
                // println!("ref_pos: {}, asm_pos: {}, score: {}", ref_pos, asm_pos, score);
                methyl_info_dict.insert((ref_pos, asm_pos), (score, node_coverage));
            }
            spos += node_seq.len();
        }
        let _ = write_methyl_bed(&methyl_info_dict, output_prefix, *hap_index, chromosome);
    }
}

pub fn start(
    graph_filename: &PathBuf,
    germline_only: bool,
    haplotype_number: usize,
    output_prefix: &PathBuf,
    het_fold_threshold: f64,
) -> AnyhowResult<()> {
    let (node_info, edge_info) = load_graph(graph_filename).unwrap();
    info!(
        "Traversing graph with germline_only: {}, hap_number: {}",
        germline_only, haplotype_number
    );

    let (haplotype_reads, node_haplotype) =
        find_node_haplotype(&node_info, haplotype_number, het_fold_threshold);
    let all_paths = enumerate_all_paths_with_haplotype(
        &node_info,
        &edge_info,
        &node_haplotype,
        haplotype_number,
    )
    .expect("Failed to enumerate all paths");
    let allseq = construct_sequences_from_haplotype_path(&node_info, &all_paths);
    // println!("allseq: {:?}", allseq.iter().map(|(index, pathlist)| format!("index: {}, path num: {}", index, pathlist.len())).collect::<Vec<_>>().join("\n"));
    let primary_haplotypes = find_full_range_haplotypes(&node_info, &node_haplotype, &allseq);
    info!("All sequences constructed: {}", primary_haplotypes.len());
    // call methylation
    call_methylation(&node_info, primary_haplotypes.clone(), output_prefix);
    info!(
        "Haplotype specific methylation signals exported: {}",
        primary_haplotypes.len()
    );
    // write assemblies
    let output_filename = PathBuf::from(format!("{}.fasta", output_prefix.to_string_lossy()));
    let _ = write_graph_path_fasta(&primary_haplotypes, &output_filename);
    info!(
        "All sequences written to fasta: {}",
        output_prefix.to_str().unwrap()
    );
    Ok(())
}
