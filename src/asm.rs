use crate::eval;
use crate::util;
use anyhow::Result as AnyhowResult;
use flate2::read::GzDecoder;
use log::{info, warn};
use ndarray::Array2;
use serde_json::Value;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::path::PathBuf;
use ndarray::s;
use rayon::prelude::*;

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
                    .unwrap_or_else(|_| panic!("Support reads not found for node: {}, {}",
                        name,
                        value["support_reads"])),
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

/// Minimum supporting reads required for a node to be considered a haplotype candidate.
const HET_MIN_NODE_SUPPORT: usize = 2;

/// Jaccard similarity |A ∩ B| / |A ∪ B|; returns 0.0 if both sets are empty.
fn read_set_jaccard(a: &HashSet<String>, b: &HashSet<String>) -> f64 {
    let union = a.union(b).count();
    if union == 0 {
        0.0
    } else {
        a.intersection(b).count() as f64 / union as f64
    }
}

/// Maximum pairwise Jaccard across all node pairs in the subset (0.0 if subset has < 2 nodes).
fn max_pairwise_jaccard(read_sets: &[HashSet<String>]) -> f64 {
    let mut worst = 0.0f64;
    for i in 0..read_sets.len() {
        for j in (i + 1)..read_sets.len() {
            worst = worst.max(read_set_jaccard(&read_sets[i], &read_sets[j]));
        }
    }
    worst
}

/// Balance score: min/max support across the subset (1.0 = perfect balance, 0.0 = degenerate).
fn min_max_support_ratio(supports: &[usize]) -> f64 {
    let max = *supports.iter().max().unwrap_or(&0) as f64;
    let min = *supports.iter().min().unwrap_or(&0) as f64;
    if max == 0.0 {
        0.0
    } else {
        min / max
    }
}

/// Identify a subset of parallel nodes at each interval that look like distinct haplotypes.
///
/// Among all feasible subsets of the largest valid size, the one maximizing the composite score
/// `(1 - max_jaccard, min_max_ratio, total_unique_reads)` is selected (lexicographic).
pub fn identify_heterozygous_nodes(
    node_info: &HashMap<String, NodeInfo>,
    hap_number: usize,
) -> HashMap<String, HashSet<String>> {
    let mut heterozygous_nodes: HashMap<String, HashSet<String>> = HashMap::new();
    if hap_number < 2 {
        return heterozygous_nodes;
    }

    let all_nodes = find_parallele_nodes(node_info);
    let mut interval_list: Vec<&String> = all_nodes.keys().collect();
    interval_list.sort_by(|a, b| {
        util::split_locus(a.to_string())
            .1
            .cmp(&util::split_locus(b.to_string()).1)
    });

    for interval_name in interval_list.iter() {
        let node_vec = all_nodes.get(*interval_name).unwrap();

        let mut candidates: Vec<(String, HashSet<String>)> = node_vec
            .iter()
            .map(|node| (node.clone(), get_read_name_list(node_info, node.clone())))
            .filter(|(_, reads)| reads.len() >= HET_MIN_NODE_SUPPORT)
            .collect();
        if candidates.len() < 2 {
            continue;
        }
        candidates.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

        let max_k = hap_number.min(candidates.len());
        let mut chosen: Option<Vec<String>> = None;

        for k in (2..=max_k).rev() {
            let mut best: Option<(f64, f64, usize, Vec<String>)> = None;

            for combo in (0..candidates.len()).combinations(k) {
                let reads: Vec<&HashSet<String>> =
                    combo.iter().map(|&i| &candidates[i].1).collect();
                let supports: Vec<usize> = reads.iter().map(|r| r.len()).collect();

                let balance = min_max_support_ratio(&supports);

                let reads_owned: Vec<HashSet<String>> =
                    reads.iter().map(|&r| r.clone()).collect();
                let max_jac = max_pairwise_jaccard(&reads_owned);

                let union_reads: HashSet<&String> =
                    reads.iter().flat_map(|r| r.iter()).collect();
                let disjointness = 1.0 - max_jac;
                let total_unique = union_reads.len();

                let candidate_score = (disjointness, balance, total_unique);
                let is_better = match &best {
                    None => true,
                    Some((d, b, t, _)) => {
                        (candidate_score.0, candidate_score.1, candidate_score.2)
                            > (*d, *b, *t)
                    }
                };
                if is_better {
                    let nodes: Vec<String> =
                        combo.iter().map(|&i| candidates[i].0.clone()).collect();
                    best = Some((disjointness, balance, total_unique, nodes));
                }
            }

            if let Some((_, _, _, nodes)) = best {
                chosen = Some(nodes);
                break;
            }
        }

        if let Some(nodes) = chosen {
            heterozygous_nodes
                .entry((*interval_name).clone())
                .or_default()
                .extend(nodes);
        }
    }
    heterozygous_nodes
}

/// Maximum number of corrections (flips) allowed per node row.
/// Bounds the DP state expansion to `sum_{i=0..=k} C(n_reads, i)` candidates per row.
const MEC_K_PER_NODE: usize = 2;
/// Soft balance penalty (corrections per 1-unit imbalance |‖A‖ − ‖B‖|).
const MEC_BALANCE_WEIGHT: f64 = 0.25;
/// Hard cap on `n_reads` for the k-cMEC DP; above this the local-search heuristic kicks in.
const MEC_KCMEC_MAX_READS: usize = 256;
/// Iteration cap for the local-search heuristic.
const MEC_LOCAL_SEARCH_MAX_ITER: usize = 100;

pub fn mec_partition_reads(
    het_matrix: &Array2<f64>,
    balance_weight: f64,
) -> (Vec<usize>, Vec<usize>, usize) {
    let (n_nodes, n_reads) = het_matrix.dim();
    if n_reads <= 1 || n_nodes == 0 {
        return ((0..n_reads).collect(), Vec::new(), 0);
    }

    if n_reads <= MEC_KCMEC_MAX_READS {
        if let Some(result) = mec_kcmec_dp(het_matrix, MEC_K_PER_NODE, balance_weight) {
            return result;
        }
    }
    mec_local_search(het_matrix, balance_weight)
}


fn mec_kcmec_dp(
    het_matrix: &Array2<f64>,
    k: usize,
    balance_weight: f64,
) -> Option<(Vec<usize>, Vec<usize>, usize)> {
    let (n_nodes, n_reads) = het_matrix.dim();
    if n_reads < 2 || n_nodes == 0 {
        return None;
    }
    let mut dp: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut seeded = false;

    for row_idx in 0..n_nodes {
        let a_indices: Vec<usize> = (0..n_reads)
            .filter(|&c| het_matrix[[row_idx, c]] > 0.5)
            .collect();
        let b_indices: Vec<usize> = (0..n_reads)
            .filter(|&c| het_matrix[[row_idx, c]] < -0.5)
            .collect();

        // Skip rows that carry no phasing information.
        if a_indices.is_empty() || b_indices.is_empty() {
            continue;
        }

        // Active = sorted union of both allele-index sets.
        let mut active: Vec<usize> =
            a_indices.iter().chain(b_indices.iter()).copied().collect();
        active.sort_unstable();

        // Enumerate all k-corrections over active reads.
        // Natural assignment: a_reads → 0 (group A), b_reads → 1 (group B), missing → 0.
        let mut local: HashMap<Vec<u8>, usize> = HashMap::new();
        for flips in 0..=k {
            for combo in active.iter().copied().combinations(flips) {
                let mut canon = vec![0u8; n_reads];
                for &c in &b_indices {
                    canon[c] = 1;
                }
                for c in &combo {
                    canon[*c] ^= 1;
                }
                // Ensure the bipartition is non-trivial over active reads.
                let active_in_b = active.iter().filter(|&&c| canon[c] == 1).count();
                if active_in_b == 0 || active_in_b == active.len() {
                    continue;
                }
                // Canonicalize: fix canon[0] = 0 (read 0 always in group A).
                if canon[0] == 1 {
                    for bit in canon.iter_mut() {
                        *bit ^= 1;
                    }
                }
                local
                    .entry(canon)
                    .and_modify(|c| {
                        if flips < *c {
                            *c = flips;
                        }
                    })
                    .or_insert(flips);
            }
        }

        if local.is_empty() {
            return None;
        }

        if !seeded {
            dp = local;
            seeded = true;
        } else {
            let mut new_dp: HashMap<Vec<u8>, usize> =
                HashMap::with_capacity(dp.len().min(local.len()));
            for (canon, prev_cost) in &dp {
                if let Some(&local_cost) = local.get(canon) {
                    new_dp.insert(canon.clone(), prev_cost + local_cost);
                }
            }
            if new_dp.is_empty() {
                return None;
            }
            dp = new_dp;
        }
    }

    if !seeded {
        // Every row lacked phasing information; return a trivial split.
        return Some(((0..n_reads).collect(), Vec::new(), 0));
    }

    let (best_canon, &best_cost) = dp.iter().min_by(|x, y| {
        let px = x.0.iter().filter(|&&bit| bit == 1).count();
        let py = y.0.iter().filter(|&&bit| bit == 1).count();
        let bal_x = balance_weight * (n_reads as f64 - 2.0 * px as f64).abs();
        let bal_y = balance_weight * (n_reads as f64 - 2.0 * py as f64).abs();
        ((*x.1 as f64) + bal_x)
            .partial_cmp(&((*y.1 as f64) + bal_y))
            .unwrap_or(std::cmp::Ordering::Equal)
    })?;

    let group_a: Vec<usize> = (0..n_reads).filter(|&i| best_canon[i] == 0).collect();
    let group_b: Vec<usize> = (0..n_reads).filter(|&i| best_canon[i] == 1).collect();
    Some((group_a, group_b, best_cost))
}

/// Compute the total number of minority-vote corrections for the current assignment,
/// counting only non-missing (non-zero) reads.  O(n_rows · n_reads).
fn mec_ternary_corrections(
    het_matrix: &Array2<f64>,
    assignment: &[u8],
) -> usize {
    let (n_nodes, n_reads) = het_matrix.dim();
    let mut corrections = 0usize;
    for r in 0..n_nodes {
        let (mut a_pos, mut a_neg, mut b_pos, mut b_neg) =
            (0usize, 0usize, 0usize, 0usize);
        for c in 0..n_reads {
            let vote = het_matrix[[r, c]];
            if vote.abs() < 0.5 {
                continue; // missing — does not contribute
            }
            let is_pos = vote > 0.0;
            match assignment[c] {
                0 => {
                    if is_pos { a_pos += 1; } else { a_neg += 1; }
                }
                _ => {
                    if is_pos { b_pos += 1; } else { b_neg += 1; }
                }
            }
        }
        corrections += a_pos.min(a_neg) + b_pos.min(b_neg);
    }
    corrections
}

/// Greedy local search: alternates between (a) majority-vote consensus per group ignoring
/// missing entries, (b) re-assigning each read to the closer haplotype (Hamming distance
/// over active votes + balance gradient).  Converges quickly but is not globally optimal.
fn mec_local_search(
    het_matrix: &Array2<f64>,
    balance_weight: f64,
) -> (Vec<usize>, Vec<usize>, usize) {
    let (n_nodes, n_reads) = het_matrix.dim();

    // Initial bipartition: alternating A/B by column index.
    let mut assignment: Vec<u8> = (0..n_reads).map(|i| (i & 1) as u8).collect();

    let mut last_objective = f64::MAX;
    for _ in 0..MEC_LOCAL_SEARCH_MAX_ITER {
        // ── Step A: consensus haplotype per group ─────────────────────────────
        // hap_X[r] = 0 means group X uses allele-A (+1) at interval r,
        //           = 1 means group X uses allele-B (-1).
        // Only non-missing (|vote| > 0.5) reads vote.
        let mut hap_a = vec![0u8; n_nodes];
        let mut hap_b = vec![0u8; n_nodes];
        for r in 0..n_nodes {
            let (mut a_pos, mut a_neg, mut b_pos, mut b_neg) =
                (0usize, 0usize, 0usize, 0usize);
            for c in 0..n_reads {
                let vote = het_matrix[[r, c]];
                if vote.abs() < 0.5 {
                    continue;
                }
                let is_pos = vote > 0.0;
                match assignment[c] {
                    0 => {
                        if is_pos { a_pos += 1; } else { a_neg += 1; }
                    }
                    _ => {
                        if is_pos { b_pos += 1; } else { b_neg += 1; }
                    }
                }
            }
            hap_a[r] = (a_neg > a_pos) as u8;
            hap_b[r] = (b_neg > b_pos) as u8;
        }

        // ── Step B: reassign each read to the closer haplotype ────────────────
        let size_a = assignment.iter().filter(|&&g| g == 0).count();
        let size_b = n_reads - size_a;
        let mut changed = false;
        for c in 0..n_reads {
            let mut dist_a = 0i64;
            let mut dist_b = 0i64;
            for r in 0..n_nodes {
                let vote = het_matrix[[r, c]];
                if vote.abs() < 0.5 {
                    continue; // missing — skip
                }
                let read_allele = if vote > 0.0 { 0u8 } else { 1u8 };
                if read_allele != hap_a[r] { dist_a += 1; }
                if read_allele != hap_b[r] { dist_b += 1; }
            }
            let current = assignment[c];
            // Balance gradient: moving from the oversize group is cheaper.
            let bal_a = (size_a as i64 - size_b as i64 - if current == 0 { 0 } else { 2 }).abs();
            let bal_b = (size_a as i64 - size_b as i64 + if current == 1 { 0 } else { 2 }).abs();
            let score_a = dist_a as f64 + balance_weight * bal_a as f64;
            let score_b = dist_b as f64 + balance_weight * bal_b as f64;
            let new_g = if score_a <= score_b { 0u8 } else { 1u8 };
            if new_g != current {
                assignment[c] = new_g;
                changed = true;
            }
        }

        // ── Step C: compute cost; bail early if converged ────────────────────
        let corrections = mec_ternary_corrections(het_matrix, &assignment);
        let sa = assignment.iter().filter(|&&g| g == 0).count();
        let sb = n_reads - sa;
        let objective = corrections as f64 + balance_weight * (sa as f64 - sb as f64).abs();

        if !changed || objective >= last_objective {
            return (
                (0..n_reads).filter(|&c| assignment[c] == 0).collect(),
                (0..n_reads).filter(|&c| assignment[c] == 1).collect(),
                corrections,
            );
        }
        last_objective = objective;
    }

    let corrections = mec_ternary_corrections(het_matrix, &assignment);
    (
        (0..n_reads).filter(|&c| assignment[c] == 0).collect(),
        (0..n_reads).filter(|&c| assignment[c] == 1).collect(),
        corrections,
    )
}

/// Build a per-interval ternary matrix for MEC phasing.
fn construct_het_interval_matrix(
    node_info: &HashMap<String, NodeInfo>,
    heterozygous_nodes: &HashMap<String, HashSet<String>>,
) -> (Array2<f64>, Vec<String>) {
    // Collect all reads across het intervals.
    let mut all_reads: HashSet<String> = HashSet::new();
    for nodes in heterozygous_nodes.values() {
        for node in nodes {
            all_reads.extend(get_read_name_list(node_info, node.clone()));
        }
    }
    let mut read_list: Vec<String> = all_reads.into_iter().collect();
    read_list.sort_unstable();

    let read_index: HashMap<&str, usize> = read_list
        .iter()
        .enumerate()
        .map(|(i, r)| (r.as_str(), i))
        .collect();

    // Sort intervals by start position for a deterministic row order.
    let mut intervals: Vec<&String> = heterozygous_nodes.keys().collect();

    intervals.sort_by(|a, b| {
        util::split_locus((*a).clone())
            .1
            .cmp(&util::split_locus((*b).clone()).1)
    });

    let n_reads = read_list.len();
    let n_rows = intervals.len();
    let mut matrix = Array2::<f64>::zeros((n_rows, n_reads));

    for (row, interval) in intervals.iter().enumerate() {
        let nodes = heterozygous_nodes.get(*interval).unwrap();
        // Sort nodes so allele-0 vs allele-1 assignment is deterministic.
        let mut sorted_nodes: Vec<&String> = nodes.iter().collect();
        sorted_nodes.sort();

        for (allele_idx, node) in sorted_nodes.iter().enumerate() {
            let vote = if allele_idx == 0 { 1.0_f64 } else { -1.0_f64 };
            for read in get_read_name_list(node_info, (*node).clone()) {
                if let Some(&col) = read_index.get(read.as_str()) {
                    if matrix[[row, col]] == 0.0 {
                        matrix[[row, col]] = vote;
                    } else {
                        // Read present in multiple alleles at this interval → conflict.
                        matrix[[row, col]] = 0.0;
                    }
                }
            }
        }
    }

    (matrix, read_list)
}

pub fn assign_haplotype_reads(
    node_info: &HashMap<String, NodeInfo>,
    heterozygous_nodes: &HashMap<String, HashSet<String>>,
    hap_number: usize,
) -> HashMap<usize, HashSet<String>> {
    let mut haplotype_reads: HashMap<usize, HashSet<String>> = HashMap::new();
    if hap_number < 2 {
        return haplotype_reads;
    }
    if !heterozygous_nodes.values().any(|s| !s.is_empty()) {
        return haplotype_reads;
    }

    // Build the ternary per-interval matrix; read_list matches the column order.
    let (het_matrix, read_list) =
        construct_het_interval_matrix(node_info, heterozygous_nodes);
    if read_list.is_empty() {
        return haplotype_reads;
    }

    let (group_a, group_b, corrections) =
        mec_partition_reads(&het_matrix, MEC_BALANCE_WEIGHT);
    info!(
        "MEC partition: |A| = {}, |B| = {}, corrections = {}",
        group_a.len(),
        group_b.len(),
        corrections
    );

    // #endregion

    haplotype_reads.insert(
        0,
        group_a.into_iter().map(|i| read_list[i].clone()).collect(),
    );
    haplotype_reads.insert(
        1,
        group_b.into_iter().map(|i| read_list[i].clone()).collect(),
    );
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

const MAX_NODE_PER_CHUNK: usize = 200;
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
    let node_list_filtered: Vec<String> = if matrix.shape()[0] > 2000 {
        let n_rows = matrix.shape()[0];
        let chunk_num = n_rows / MAX_NODE_PER_CHUNK;
        let chunks_filtered: Vec<String> = (0..chunk_num)
            .into_par_iter()
            .flat_map_iter(|i| {
                let raw_start = i * MAX_NODE_PER_CHUNK;
                let (start, end) = if raw_start + MAX_NODE_PER_CHUNK > n_rows {
                    (n_rows - MAX_NODE_PER_CHUNK, n_rows)
                } else {
                    (raw_start, raw_start + MAX_NODE_PER_CHUNK)
                };
                let chunk = matrix.slice(s![start..end, ..]).to_owned();
                let node_list_chunk = node_list[start..end].to_vec();
                util::permutation_test(&chunk, 0.5, 100, node_list_chunk).into_iter()
            })
            .collect();
        chunks_filtered.into_iter().collect::<HashSet<_>>().into_iter().collect()
    } else {
        util::permutation_test(&matrix, 0.5, 100, node_list.clone())
    };
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
    haplotype_reads: &HashMap<usize, HashSet<String>>
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
        let mut haplotype_intersection_src = node_haplotype
            .get(&src)
            .unwrap_or(&HashSet::new())
            .intersection(&haplotype_index)
            .cloned()
            .collect::<HashSet<_>>();
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
    // println!("current_node: {:?}, haplotype_index: {:?}", current_node, haplotype_index);
    let next_nodes = edge_info.get(current_node).unwrap();
    for next_node in next_nodes {
        let haplotypes_next = if node_haplotype.contains_key(next_node) {
            node_haplotype.get(&next_node.clone()).unwrap()
        } else {
            &HashSet::new()
        };
        let mut haplotype_intersection_clone = haplotype_index
            .intersection(haplotypes_next)
            .cloned()
            .collect::<HashSet<_>>();
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
            let node_info_dict = node_info.get(node).unwrap_or_else(|| panic!("Node {} not found in node_info, path Number: {:?}",
                node, path));
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
            all_sequences
                .entry(*hap_ind)
                .or_insert(Vec::new())
                .push((path.clone(), sequence.clone(), read_names.clone()));
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
    all_sequences: &HashMap<usize, Vec<(Vec<String>, String, HashSet<String>)>>,
) -> HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)> {
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
        // first compare the 4th element, then the 3th element
        // full_sequences.sort_by(|a, b| b.3.cmp(&a.3).then(b.4.cmp(&a.4)));
        full_sequences.sort_by(|a, b| b.4.cmp(&a.4).then(b.3.cmp(&a.3)));
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

    for (hap, node_list) in haplotype_nodes.iter() {
        let interval_node = find_parallele_nodes_from_nodelist(
            &node_list.iter().cloned().collect::<Vec<_>>(),
        );
        let mut interval_list = interval_node.keys().collect::<Vec<_>>();
        interval_list.sort(); // ascending order a < b

        let read_list = haplotype_reads.get(hap).unwrap().clone();
        // #region agent log (H3: filter_haplotype_nodes tie-break)

        // #endregion
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
                filtered_haplotype_nodes
                    .entry(*hap)
                    .or_insert(HashSet::new())
                    .insert(best_node.clone());
            } else {
                warn!(
                    "hap: {}, interval: {}, nodes: {:?}",
                    hap,
                    interval,
                    nodes
                        .clone()
                        .iter().cloned()
                        .collect::<Vec<_>>()
                        .join(", ")
                );
            }
        }

        // #endregion
    }
    filtered_haplotype_nodes
}

pub fn find_node_haplotype(
    node_info: &HashMap<String, NodeInfo>,
    hap_number: usize,
) -> (
    HashMap<usize, HashSet<String>>,
    HashMap<String, HashSet<usize>>,
) {
    if hap_number == 1 {
        let node_haplotype = find_most_supported_path(node_info);
        return (HashMap::new(), node_haplotype);
    }else if hap_number == 2 {
        let heterozygous_nodes = identify_heterozygous_nodes(node_info, hap_number);
        info!("heterozygous_nodes: {:?}", heterozygous_nodes.len());
        // let filtered_heterozygous_nodes = filter_heterozygous_nodes(node_info, &heterozygous_nodes);
        // info!(
        //     "filtered_heterozygous_nodes: {:?}",
        //     filtered_heterozygous_nodes.len()
        // );

        let haplotype_reads =
            assign_haplotype_reads(node_info, &heterozygous_nodes, hap_number);

        // let (haplotype_reads_new, haplotype_nodes_new) =
        //     assign_unassigned_reads(node_info, &haplotype_reads);
        // info!(
        //     "haplotype_reads: {:?}",
        //     haplotype_reads
        //         .iter()
        //         .map(|(hap, reads)| format!("hap: {}, reads: {}", hap, reads.len()))
        //         .collect::<Vec<_>>()
        //         .join(", ")
        // );
        // info!(
        //     "haplotype_reads_new: {:?}",
        //     haplotype_reads_new
        //         .iter()
        //         .map(|(hap, reads)| format!("hap: {}, reads: {}", hap, reads.len()))
        //         .collect::<Vec<_>>()
        //         .join(", ")
        // );
        // println!("haplotype_reads_new: {}, {}, {:?}", haplotype_reads_new.get(&0).unwrap().len(), haplotype_reads_new.get(&1).unwrap().len(), haplotype_reads_new.get(&0).unwrap().intersection(haplotype_reads_new.get(&1).unwrap()).count());
        let mut total_reads = HashSet::new();
        for node in node_info.keys(){
            let read_names = get_read_name_list(node_info, node.clone());
            total_reads.extend(read_names);
        }
        info!(
            "total reads: {:?}",
            total_reads.len()
        );        
        if haplotype_reads.is_empty() {
            let node_haplotype = find_most_supported_path(node_info);
            return (haplotype_reads, node_haplotype)
        } else {
            let mut haplotype_nodes = HashMap::new();
            let read_to_nodes = assign_node_to_reads(node_info);
            
            for (hap, reads) in haplotype_reads.iter() {
                for read in reads.iter() {
                    let nodes_on_read = read_to_nodes.get(read).unwrap().clone();
                    haplotype_nodes.entry(*hap).or_insert(HashSet::new()).extend(nodes_on_read);
                }
            }
            let filtered_haplotype_nodes =
                filter_haplotype_nodes(node_info, &haplotype_nodes, &haplotype_reads);

            let node_haplotype = assign_haplotype_to_nodes(&filtered_haplotype_nodes);
            
            // #region agent log (H4+H5: node_haplotype label distribution)
            {
                let mut nodes_hap0 = 0usize;
                let mut nodes_hap1 = 0usize;
                let mut nodes_both = 0usize;
                for haps in node_haplotype.values() {
                    let has0 = haps.contains(&0);
                    let has1 = haps.contains(&1);
                    if has0 && has1 { nodes_both += 1; }
                    else if has0 { nodes_hap0 += 1; }
                    else if has1 { nodes_hap1 += 1; }
                }
                let reads_hap0 = haplotype_reads.get(&0).map(|s| s.len()).unwrap_or(0);
                let reads_hap1 = haplotype_reads.get(&1).map(|s| s.len()).unwrap_or(0);
            }
            // #endregion
        return (haplotype_reads, node_haplotype)
    }}else{
        warn!("Unsupported hap_number: {}", hap_number);
        (HashMap::new(), HashMap::new())
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
            let node_seq = node_info.get(node).unwrap().seq.clone();
            let node_coverage = node_info.get(node).unwrap().support_reads;
            for (pos, score) in methyl_info.into_iter() {
                let asm_pos = pos + spos;
                let ref_pos = *position_mapping.get(&pos).unwrap_or(&0);
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
    output_prefix: &PathBuf
) -> AnyhowResult<(HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)>, HashMap<String, NodeInfo>, HashMap<String, Vec<String>>)> {
    let (node_info, edge_info) = load_graph(graph_filename).unwrap();
    info!(
        "Traversing graph with germline_only: {}, hap_number: {}",
        germline_only, haplotype_number
    );

    let (haplotype_reads, node_haplotype) =
        find_node_haplotype(&node_info, haplotype_number);
        
    // println!("node_haplotype: {:?}", node_haplotype);
    let all_paths = enumerate_all_paths_with_haplotype(
        &node_info,
        &edge_info,
        &node_haplotype,
        haplotype_number,
    )
    .expect("Failed to enumerate all paths");
    let allseq = construct_sequences_from_haplotype_path(&node_info, &all_paths);
    let primary_haplotypes =
        find_full_range_haplotypes(&node_info, &allseq);
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
    Ok((primary_haplotypes, node_info.clone(), edge_info.clone()))
}

#[cfg(test)]
mod identify_heterozygous_tests {
    use super::*;

    fn make_node(reads: &[&str]) -> NodeInfo {
        NodeInfo {
            seq: "A".to_string(),
            cigar: "10=".to_string(),
            support_reads: reads.len(),
            allele_frequency: "0.5".to_string(),
            read_names: reads.join(","),
            methyl_info: HashMap::new(),
        }
    }

    fn insert(map: &mut HashMap<String, NodeInfo>, id: &str, reads: &[&str]) {
        map.insert(id.to_string(), make_node(reads));
    }

    #[test]
    fn balanced_disjoint_pair_accepted() {
        let mut nodes = HashMap::new();
        insert(&mut nodes, "H.chr1:1-100.0", &["r1", "r2", "r3", "r4"]);
        insert(&mut nodes, "H.chr1:1-100.1", &["r5", "r6", "r7", "r8"]);

        let het = identify_heterozygous_nodes(&nodes, 2);
        let chosen = het.get("chr1:1-100").expect("interval should be heterozygous");
        assert_eq!(chosen.len(), 2);
    }

    #[test]
    fn overlapping_reads_rejected_for_full_haps() {
        let mut nodes = HashMap::new();
        insert(&mut nodes, "H.chr1:1-100.0", &["r1", "r2", "r3", "r4"]);
        insert(&mut nodes, "H.chr1:1-100.1", &["r1", "r2", "r3", "r5"]);

        let het = identify_heterozygous_nodes(&nodes, 2);
        assert!(het.get("chr1:1-100").is_none());
    }

    #[test]
    fn imbalanced_support_rejected() {
        let mut nodes = HashMap::new();
        insert(&mut nodes, "H.chr1:1-100.0", &["r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8"]);
        insert(&mut nodes, "H.chr1:1-100.1", &["r9", "r10"]);

        let het = identify_heterozygous_nodes(&nodes, 2);
        assert!(het.get("chr1:1-100").is_none());
    }

    #[test]
    fn low_support_node_skipped() {
        let mut nodes = HashMap::new();
        insert(&mut nodes, "H.chr1:1-100.0", &["r1", "r2", "r3"]);
        insert(&mut nodes, "H.chr1:1-100.1", &["r4"]);

        let het = identify_heterozygous_nodes(&nodes, 2);
        assert!(het.get("chr1:1-100").is_none());
    }

    #[test]
    fn back_off_from_three_to_two() {
        let mut nodes = HashMap::new();
        insert(&mut nodes, "H.chr1:1-100.0", &["r1", "r2", "r3", "r4"]);
        insert(&mut nodes, "H.chr1:1-100.1", &["r5", "r6", "r7", "r8"]);
        // Third allele shares heavily with both other alleles → no valid 3-subset.
        insert(&mut nodes, "H.chr1:1-100.2", &["r1", "r2", "r5", "r6"]);

        let het = identify_heterozygous_nodes(&nodes, 3);
        let chosen = het.get("chr1:1-100").expect("partial het set should be returned");
        assert_eq!(chosen.len(), 2, "should back off to disjoint pair");
        assert!(chosen.contains("H.chr1:1-100.0"));
        assert!(chosen.contains("H.chr1:1-100.1"));
    }

    #[test]
    fn picks_disjoint_subset_among_many_candidates() {
        let mut nodes = HashMap::new();
        // Two pairs are roughly balanced; only one pair is also disjoint.
        insert(&mut nodes, "H.chr1:1-100.0", &["r1", "r2", "r3", "r4"]);
        insert(&mut nodes, "H.chr1:1-100.1", &["r1", "r2", "r3", "r4"]);
        insert(&mut nodes, "H.chr1:1-100.2", &["r5", "r6", "r7", "r8"]);

        let het = identify_heterozygous_nodes(&nodes, 2);
        let chosen = het.get("chr1:1-100").expect("should find disjoint pair");
        assert!(chosen.contains("H.chr1:1-100.2"));
        let other = chosen
            .iter()
            .find(|id| id.as_str() != "H.chr1:1-100.2")
            .expect("expected a second allele");
        assert!(other == "H.chr1:1-100.0" || other == "H.chr1:1-100.1");
    }

    #[test]
    fn mec_recovers_clean_diploid_partition() {
        // Ternary matrix: +1 = votes allele-A, -1 = votes allele-B, 0 = missing.
        // Two intervals, 3 reads per haplotype, error-free → 0 corrections expected.
        let m = ndarray::array![
            // reads:   0     1     2     3     4     5
            [ 1.0,  1.0,  1.0, -1.0, -1.0, -1.0], // interval 0
            [ 1.0,  1.0,  1.0, -1.0, -1.0, -1.0], // interval 1
        ];
        let (a, b, corr) = mec_partition_reads(&m, 0.0);
        assert_eq!(corr, 0);
        let mut a_sorted = a.clone();
        a_sorted.sort();
        let mut b_sorted = b.clone();
        b_sorted.sort();
        assert!(
            (a_sorted == vec![0, 1, 2] && b_sorted == vec![3, 4, 5])
                || (a_sorted == vec![3, 4, 5] && b_sorted == vec![0, 1, 2])
        );
    }

    #[test]
    fn mec_corrects_single_flip() {
        // Read 0 votes allele-B at interval 0 but allele-A at interval 1.
        // One correction needed to assign read 0 to group A.
        let m = ndarray::array![
            [-1.0,  1.0,  1.0, -1.0, -1.0, -1.0], // interval 0: read 0 flipped
            [ 1.0,  1.0,  1.0, -1.0, -1.0, -1.0], // interval 1: clean
        ];
        let (_, _, corr) = mec_partition_reads(&m, 0.0);
        assert_eq!(corr, 1, "exactly one entry should need flipping");
    }

    #[test]
    fn mec_balance_penalty_breaks_ties() {
        // 1 interval, reads 0-1 vote allele-A, reads 2-3 vote allele-B.
        // Zero-correction 2-2 split should beat any 3-1 split with balance_weight = 1.
        let m = ndarray::array![[1.0, 1.0, -1.0, -1.0]];
        let (a, b, corr) = mec_partition_reads(&m, 1.0);
        assert_eq!(corr, 0);
        assert_eq!(a.len() + b.len(), 4);
        assert_eq!((a.len() as i64 - b.len() as i64).abs(), 0);
    }

    #[test]
    fn mec_missing_reads_are_wildcards() {
        // Interval 0: reads 0-2 vote allele-A, reads 3-4 vote allele-B, read 5 missing.
        // Interval 1: reads 0-2 vote allele-A, reads 3-4 vote allele-B, read 5 missing.
        // Read 5 has no information → should not inflate the correction count.
        let m = ndarray::array![
            [ 1.0,  1.0,  1.0, -1.0, -1.0,  0.0],
            [ 1.0,  1.0,  1.0, -1.0, -1.0,  0.0],
        ];
        let (a, b, corr) = mec_partition_reads(&m, 0.0);
        assert_eq!(corr, 0, "missing reads must not add corrections");
        // Reads 0-2 and 3-4 should land in opposite groups; read 5 can be anywhere.
        let a_set: std::collections::HashSet<usize> = a.into_iter().collect();
        let b_set: std::collections::HashSet<usize> = b.into_iter().collect();
        let active_a: Vec<usize> = [0, 1, 2].iter().filter(|&&r| a_set.contains(&r)).copied().collect();
        let active_b: Vec<usize> = [0, 1, 2].iter().filter(|&&r| b_set.contains(&r)).copied().collect();
        assert!(
            active_a.len() == 3 || active_b.len() == 3,
            "reads 0-2 should be in the same group"
        );
    }

    #[test]
    fn mec_kcmec_falls_back_when_too_many_flips_needed() {
        // Interval 0 alternates (+/-), intervals 1-3 are block (+/+/+/-/-/-).
        // The optimal split {0,1,2}|{3,4,5} satisfies intervals 1-3 perfectly
        // and costs 2 corrections at interval 0 (reads 1 and 4).
        let m = ndarray::array![
            [ 1.0, -1.0,  1.0, -1.0,  1.0, -1.0], // interval 0: alternating
            [ 1.0,  1.0,  1.0, -1.0, -1.0, -1.0], // interval 1: block
            [ 1.0,  1.0,  1.0, -1.0, -1.0, -1.0], // interval 2
            [ 1.0,  1.0,  1.0, -1.0, -1.0, -1.0], // interval 3
        ];
        let (a, b, corr) = mec_partition_reads(&m, 0.0);
        assert_eq!(a.len() + b.len(), 6);
        assert!(corr <= 3, "should produce a reasonable partition");
    }

}