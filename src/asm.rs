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
use csv::Writer;

#[derive(Debug, Clone)]
pub struct NodeInfo {
    pub seq: String,
    pub cigar: String,
    pub support_reads: usize,
    pub allele_frequency: String,
    pub read_names: String,
    pub methyl_info: HashMap<usize, f32>,
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
        let (chromo, start, end )= util::split_locus(interval_name.to_string());
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
                .extend(nodes.clone());
        }
    }
    heterozygous_nodes
}

/// Soft balance penalty (corrections per 1-unit imbalance between cluster sizes).
const CLUSTER_BALANCE_WEIGHT: f64 = 0.25;
/// Number of farthest-first seeded restarts for the read-clustering phaser.
const CLUSTER_RESTARTS: usize = 64;
/// Iteration cap per restart for the Lloyd-style centroid refinement.
const CLUSTER_MAX_ITER: usize = 100;
/// Minimum number of co-covered het sites required to trust a read↔centroid comparison.
const CLUSTER_MIN_OVERLAP: usize = 1;

/// Sparse ternary representation of a single read: `(row index, sign ∈ {-1, +1})`,
/// kept sorted by row so two reads can be compared with a linear merge.
type ReadVec = Vec<(usize, i8)>;

/// Build per-read sparse vectors from the ternary het matrix (matrix columns = reads).
fn build_read_vectors(het_matrix: &Array2<f64>) -> Vec<ReadVec> {
    let (n_nodes, n_reads) = het_matrix.dim();
    let mut reads: Vec<ReadVec> = vec![Vec::new(); n_reads];
    for r in 0..n_nodes {
        for c in 0..n_reads {
            let v = het_matrix[[r, c]];
            if v > 0.5 {
                reads[c].push((r, 1i8));
            } else if v < -0.5 {
                reads[c].push((r, -1i8));
            }
        }
    }
    reads
}

/// Disagreement between a read and a dense ternary centroid, restricted to positions
/// where both are informative. Returns `(disagreements, overlap)`.
fn read_centroid_disagreement(read: &ReadVec, centroid: &[i8]) -> (usize, usize) {
    let mut disagree = 0usize;
    let mut overlap = 0usize;
    for &(row, sign) in read.iter() {
        let c = centroid[row];
        if c == 0 {
            continue;
        }
        overlap += 1;
        if c != sign {
            disagree += 1;
        }
    }
    (disagree, overlap)
}

/// Disagreement fraction between two reads over their co-covered het sites.
/// Returns `(fraction, overlap)`; the fraction is 1.0 when there is no overlap.
fn read_pair_disagreement(a: &ReadVec, b: &ReadVec) -> (f64, usize) {
    let (mut i, mut j) = (0usize, 0usize);
    let mut overlap = 0usize;
    let mut disagree = 0usize;
    while i < a.len() && j < b.len() {
        let (ra, sa) = a[i];
        let (rb, sb) = b[j];
        match ra.cmp(&rb) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                overlap += 1;
                if sa != sb {
                    disagree += 1;
                }
                i += 1;
                j += 1;
            }
        }
    }
    if overlap == 0 {
        (1.0, 0)
    } else {
        (disagree as f64 / overlap as f64, overlap)
    }
}

/// Recompute a dense ternary centroid for `members` by per-row majority vote.
fn compute_centroid(reads: &[ReadVec], members: &[usize], n_nodes: usize) -> Vec<i8> {
    let mut score = vec![0i32; n_nodes];
    for &m in members {
        for &(row, sign) in reads[m].iter() {
            score[row] += sign as i32;
        }
    }
    score
        .into_iter()
        .map(|s| match s.cmp(&0) {
            std::cmp::Ordering::Greater => 1i8,
            std::cmp::Ordering::Less => -1i8,
            std::cmp::Ordering::Equal => 0i8,
        })
        .collect()
}

/// Total minimum-error corrections for a clustering: for every row and every cluster, the
/// minority-vote count, summed over all rows and clusters.
fn cluster_corrections(
    reads: &[ReadVec],
    assignment: &[usize],
    n_clusters: usize,
    n_nodes: usize,
) -> usize {
    let mut pos = vec![vec![0u32; n_nodes]; n_clusters];
    let mut neg = vec![vec![0u32; n_nodes]; n_clusters];
    for (read_idx, &cl) in assignment.iter().enumerate() {
        for &(row, sign) in reads[read_idx].iter() {
            if sign > 0 {
                pos[cl][row] += 1;
            } else {
                neg[cl][row] += 1;
            }
        }
    }
    let mut corrections = 0usize;
    for cl in 0..n_clusters {
        for row in 0..n_nodes {
            corrections += pos[cl][row].min(neg[cl][row]) as usize;
        }
    }
    corrections
}

/// Farthest-first seeding: pick `hap_number` reads as initial centroids, starting from
/// `first` and repeatedly adding the read whose pattern is most dissimilar (largest minimum
/// disagreement) to the reads already chosen. Different `first` values let restarts explore
/// distinct basins.
fn seed_clusters(reads: &[ReadVec], hap_number: usize, first: usize) -> Vec<usize> {
    let n_reads = reads.len();
    let mut seeds = vec![first];
    while seeds.len() < hap_number {
        let mut best_read: Option<usize> = None;
        let mut best_key = (f64::MIN, 0usize);
        for c in 0..n_reads {
            if reads[c].is_empty() || seeds.contains(&c) {
                continue;
            }
            let mut min_frac = f64::MAX;
            let mut acc_overlap = 0usize;
            let mut comparable = false;
            for &s in seeds.iter() {
                let (frac, ov) = read_pair_disagreement(&reads[c], &reads[s]);
                if ov >= CLUSTER_MIN_OVERLAP {
                    comparable = true;
                    acc_overlap += ov;
                    if frac < min_frac {
                        min_frac = frac;
                    }
                }
            }
            if !comparable {
                continue;
            }
            // Prefer reads far from existing seeds, breaking ties by larger overlap.
            let key = (min_frac, acc_overlap);
            if key > best_key {
                best_key = key;
                best_read = Some(c);
            }
        }
        match best_read {
            Some(c) => seeds.push(c),
            None => break, // not enough informative reads to seed another cluster
        }
    }
    seeds
}

/// Partition reads into `hap_number` haplotype clusters from the ternary het matrix.
///
/// Reads are clustered by the similarity of their `{-1, 0, +1}` patterns: each read is
/// compared against a per-cluster consensus spanning the *whole* interval, which makes the
/// assignment robust to switch errors (the bipartition can no longer silently flip
/// orientation midway, as happened with the row-by-row DP). Farthest-first seeding plus many
/// restarts, scored by total minimum-error corrections, avoids degenerate solutions where
/// both haplotypes collapse onto a single truth haplotype.
///
/// Returns the clusters (read column indices, largest first) and the total MEC corrections.
pub fn cluster_partition_reads(
    het_matrix: &Array2<f64>,
    hap_number: usize,
    balance_weight: f64,
) -> (Vec<Vec<usize>>, usize) {
    let (n_nodes, n_reads) = het_matrix.dim();
    let k = hap_number.max(1);

    let trivial = || {
        let mut clusters = vec![Vec::new(); k];
        clusters[0] = (0..n_reads).collect::<Vec<_>>();
        (clusters, 0usize)
    };
    if n_reads <= 1 || n_nodes == 0 || k == 1 {
        return trivial();
    }

    let reads = build_read_vectors(het_matrix);

    // Seed restarts from the highest-coverage reads first (most informative sites).
    let mut informative: Vec<usize> = (0..n_reads).filter(|&c| !reads[c].is_empty()).collect();
    informative.sort_by(|&a, &b| reads[b].len().cmp(&reads[a].len()).then(a.cmp(&b)));
    if informative.is_empty() {
        return trivial();
    }

    let restarts = CLUSTER_RESTARTS.min(informative.len());
    let mut best: Option<(f64, usize, usize, Vec<usize>)> = None; // (objective, corrections, eff_k, assignment)

    for restart in 0..restarts {
        let first = informative[restart];
        let seeds = seed_clusters(&reads, k, first);
        if seeds.len() < 2 {
            continue;
        }
        let eff_k = seeds.len();

        let mut centroids: Vec<Vec<i8>> = seeds
            .iter()
            .map(|&s| {
                let mut dense = vec![0i8; n_nodes];
                for &(row, sign) in reads[s].iter() {
                    dense[row] = sign;
                }
                dense
            })
            .collect();

        let mut assignment = vec![0usize; n_reads];
        let mut prev_assignment: Vec<usize> = Vec::new();

        for _iter in 0..CLUSTER_MAX_ITER {
            let mut sizes = vec![0usize; eff_k];

            // Assign informative reads to the closest consensus over the whole interval.
            for &c in informative.iter() {
                let mut best_cl = 0usize;
                let mut best_frac = f64::MAX;
                let mut best_overlap = 0usize;
                let mut found = false;
                for cl in 0..eff_k {
                    let (dis, ov) = read_centroid_disagreement(&reads[c], &centroids[cl]);
                    if ov < CLUSTER_MIN_OVERLAP {
                        continue;
                    }
                    let frac = dis as f64 / ov as f64;
                    let better = !found
                        || frac < best_frac - 1e-9
                        || ((frac - best_frac).abs() <= 1e-9 && ov > best_overlap);
                    if better {
                        found = true;
                        best_frac = frac;
                        best_overlap = ov;
                        best_cl = cl;
                    }
                }
                if !found {
                    // No overlap with any centroid this round: keep clusters balanced.
                    best_cl = (0..eff_k).min_by_key(|&cl| sizes[cl]).unwrap();
                }
                assignment[c] = best_cl;
                sizes[best_cl] += 1;
            }
            // Reads with no informative site carry no phasing signal → balance filler.
            for c in 0..n_reads {
                if reads[c].is_empty() {
                    let cl = (0..eff_k).min_by_key(|&x| sizes[x]).unwrap();
                    assignment[c] = cl;
                    sizes[cl] += 1;
                }
            }

            if assignment == prev_assignment {
                break;
            }
            prev_assignment = assignment.clone();

            let mut members: Vec<Vec<usize>> = vec![Vec::new(); eff_k];
            for (idx, &cl) in assignment.iter().enumerate() {
                members[cl].push(idx);
            }
            for cl in 0..eff_k {
                centroids[cl] = compute_centroid(&reads, &members[cl], n_nodes);
            }
        }

        let corrections = cluster_corrections(&reads, &assignment, eff_k, n_nodes);
        let mut sizes = vec![0usize; eff_k];
        for &cl in assignment.iter() {
            sizes[cl] += 1;
        }
        let imbalance = (sizes.iter().max().unwrap() - sizes.iter().min().unwrap()) as f64;
        let objective = corrections as f64 + balance_weight * imbalance;

        if best.as_ref().map_or(true, |(bo, _, _, _)| objective < *bo) {
            best = Some((objective, corrections, eff_k, assignment.clone()));
        }
    }

    let (_, corrections, _eff_k, assignment) = match best {
        Some(b) => b,
        None => return trivial(),
    };

    let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); k];
    for (idx, &cl) in assignment.iter().enumerate() {
        clusters[cl.min(k - 1)].push(idx);
    }
    // Largest cluster first for deterministic, stable haplotype indexing.
    clusters.sort_by(|a, b| b.len().cmp(&a.len()));
    (clusters, corrections)
}

fn write_matrix_to_csv<P: AsRef<Path>>(
    matrix: &Array2<f64>,
    groups: &[Vec<usize>],
    nodelist: &[String],
    read_set: &[String],
    path: P,
) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = Writer::from_writer(file);

    // Column order groups reads by their assigned cluster so the CSV visually separates haplotypes.
    let ordered_cols: Vec<usize> = groups.iter().flat_map(|g| g.iter().copied()).collect();

    let mut header = vec!["node".to_string()];
    for &col in ordered_cols.iter() {
        header.push(read_set[col].clone());
    }
    writer.write_record(&header)?;

    for (row_idx, node) in nodelist.iter().enumerate() {
        let mut row = vec![node.clone()];
        for &col in ordered_cols.iter() {
            row.push(matrix[[row_idx, col]].to_string());
        }
        writer.write_record(&row)?;
    }

    writer.flush()?;
    Ok(())
}

/// Build a per-interval ternary matrix for MEC phasing.
fn construct_het_interval_matrix(
    node_info: &HashMap<String, NodeInfo>,
    heterozygous_nodes: &HashMap<String, HashSet<String>>,
) -> (Array2<f64>, Vec<String>, Vec<String>) {
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
    let mut het_nodes: Vec<String> = Vec::new();
    for (interval, nodes) in heterozygous_nodes.iter() {
        het_nodes.extend(nodes.iter().cloned());
    }
    info!("het_nodes: {:?}", het_nodes.len());
    het_nodes.sort_by(|a, b| a.cmp(b));

    let n_reads = read_list.len();
    let n_rows = het_nodes.len();
    let mut matrix = Array2::<f64>::zeros((n_rows, n_reads));

    for (row, node) in het_nodes.iter().enumerate() {
        // let vote = if allele_idx == 0 { 1.0_f64 } else { -1.0_f64 };
        // get all read in the interval
        let interval = node.split(".").collect::<Vec<_>>()[1];
        let interval_nodes = heterozygous_nodes.get(interval).unwrap().iter().cloned().collect::<Vec<_>>();
        let mut total_reads = HashSet::new();
        for node in interval_nodes.iter(){
            let read_names = get_read_name_list(node_info, node.clone());
            total_reads.extend(read_names);
        }

        let current_read_list = get_read_name_list(node_info, node.clone());
        for read in  read_list.iter(){
            let vote = if total_reads.contains(read){
                if current_read_list.contains(read){
                    1.0
                } else {
                    -1.0
                }
            } else {
                0.0
            };
            
            if let Some(&col) = read_index.get(read.as_str()) {
                matrix[[row, col]] = vote;
            }
        }
        
    }

    (matrix, read_list, het_nodes)
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
    let (het_matrix, read_list, het_nodes) =
        construct_het_interval_matrix(node_info, heterozygous_nodes);
    if read_list.is_empty() {
        return haplotype_reads;
    }
    
    let (clusters, corrections) =
        cluster_partition_reads(&het_matrix, hap_number, CLUSTER_BALANCE_WEIGHT);
    write_matrix_to_csv(&het_matrix, &clusters, &het_nodes, &read_list, "het_matrix.csv").unwrap();

    info!(
        "Cluster partition: cluster sizes = {:?}, corrections = {}",
        clusters.iter().map(|c| c.len()).collect::<Vec<_>>(),
        corrections
    );
    for (hap, cluster) in clusters.into_iter().enumerate() {
        haplotype_reads.insert(
            hap,
            cluster.into_iter().map(|i| read_list[i].clone()).collect(),
        );
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
        current_path.pop();
    }
}

pub fn construct_sequences_from_haplotype_path(
    node_info: &HashMap<String, NodeInfo>,
    all_paths: &Vec<(Vec<String>, HashSet<usize>)>,
) -> HashMap<usize, Vec<(Vec<String>, String, HashSet<String>)>> {
    let mut all_sequences = HashMap::new();
    for (path_index, (path, haplotype_index)) in all_paths.iter().enumerate() {
        let mut sequence = String::new();
        let mut read_names = HashSet::new();

        for node in path.iter() {
            let node_info_dict = node_info.get(node).unwrap_or_else(|| panic!("Node {} not found in node_info, path Number: {:?}",
                node, path));
            let read_names_list_clone = get_read_name_list(node_info, node.clone());
            read_names.extend(read_names_list_clone.clone());
            sequence += &node_info_dict.seq;
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

        let seq_len = sequence.len();
        let full_lines = seq_len / chars_per_line;
        for i in 0..full_lines {
            let start = i * chars_per_line;
            let end = start + chars_per_line;
            writeln!(file, "{}", &sequence[start..end])?;
        }
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

#[derive(Clone)]
struct PathState {
    path: Vec<String>,
    sequence: String,
    read_names: HashSet<String>,
    supports: usize,
    span: usize,
}

fn path_state_from_node(node_info: &HashMap<String, NodeInfo>, node: &str) -> Option<PathState> {
    let info = node_info.get(node)?;
    let read_names = get_read_name_list(node_info, node.to_string());
    let (start, end) =
        eval::find_alignment_intervals([node.to_string()].iter().map(|x| x.as_str()).collect())
            .ok()?;
    Some(PathState {
        path: vec![node.to_string()],
        sequence: info.seq.clone(),
        read_names,
        supports: info.support_reads,
        span: end.saturating_sub(start),
    })
}

fn extend_path_state(
    node_info: &HashMap<String, NodeInfo>,
    prev: &PathState,
    node: &str,
) -> Option<PathState> {
    let info = node_info.get(node)?;
    let mut path = prev.path.clone();
    path.push(node.to_string());
    let mut read_names = prev.read_names.clone();
    read_names.extend(get_read_name_list(node_info, node.to_string()));
    let sequence = format!("{}{}", prev.sequence, info.seq);
    let (start, end) = eval::find_alignment_intervals(path.iter().map(|x| x.as_str()).collect()).ok()?;
    Some(PathState {
        path: path.clone(),
        sequence,
        read_names,
        supports: get_supports(node_info, &path),
        span: end.saturating_sub(start),
    })
}

fn path_state_better(a: &PathState, b: &PathState) -> bool {
    a.span > b.span || (a.span == b.span && a.supports > b.supports)
}

/// Layer-wise DP over the window graph: one best path per haplotype (max span, then supports).
pub fn find_best_haplotype_paths_dp(
    node_info: &HashMap<String, NodeInfo>,
    edge_info: &HashMap<String, Vec<String>>,
    node_haplotype: &HashMap<String, HashSet<usize>>,
    haplotype_number: usize,
) -> HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)> {
    let mut pred_map: HashMap<String, Vec<String>> = HashMap::new();
    for (src, dsts) in edge_info.iter() {
        for dst in dsts {
            pred_map.entry(dst.clone()).or_default().push(src.clone());
        }
    }

    let parallel_nodes = find_parallele_nodes(node_info);
    let mut interval_list: Vec<String> = parallel_nodes.keys().cloned().collect();
    interval_list.sort_by(|a, b| {
        util::split_locus(a.clone())
            .1
            .cmp(&util::split_locus(b.clone()).1)
    });

    let mut node_interval: HashMap<String, String> = HashMap::new();
    for (interval, nodes) in &parallel_nodes {
        for node in nodes {
            node_interval.insert(node.clone(), interval.clone());
        }
    }

    let mut best_by_hap: HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)> =
        HashMap::new();

    for hap in 0..haplotype_number {
        let mut dp: HashMap<String, PathState> = HashMap::new();

        for (interval_idx, interval) in interval_list.iter().enumerate() {
            let Some(nodes) = parallel_nodes.get(interval) else {
                continue;
            };
            let prev_interval = interval_idx.checked_sub(1).map(|i| &interval_list[i]);
            for node in nodes {
                let Some(tags) = node_haplotype.get(node) else {
                    continue;
                };
                if !tags.contains(&hap) {
                    continue;
                }

                let preds = pred_map.get(node).map(|v| v.as_slice()).unwrap_or(&[]);
                let mut best: Option<PathState> = None;

                for pred in preds {
                    if let Some(prev_int) = prev_interval {
                        if node_interval.get(pred) != Some(prev_int) {
                            continue;
                        }
                    }
                    if let Some(prev) = dp.get(pred) {
                        if let Some(extended) = extend_path_state(node_info, prev, node) {
                            if best.as_ref().map_or(true, |b| path_state_better(&extended, b)) {
                                best = Some(extended);
                            }
                        }
                    }
                }

                if best.is_none() {
                    if let Some(start) = path_state_from_node(node_info, node) {
                        best = Some(start);
                    }
                }

                if let Some(state) = best {
                    dp.entry(node.clone())
                        .and_modify(|existing| {
                            if path_state_better(&state, existing) {
                                *existing = state.clone();
                            }
                        })
                        .or_insert(state);
                }
            }
        }

        let mut hap_best: Option<PathState> = None;
        for (node, state) in &dp {
            if !node_haplotype
                .get(node)
                .map(|t| t.contains(&hap))
                .unwrap_or(false)
            {
                continue;
            }
            if hap_best.as_ref().map_or(true, |b| path_state_better(state, b)) {
                hap_best = Some(state.clone());
            }
        }

        if let Some(state) = hap_best {
            best_by_hap.insert(
                hap,
                (
                    state.path,
                    state.sequence,
                    state.read_names,
                    state.supports,
                    state.span,
                ),
            );
        }
    }

    best_by_hap
}

pub fn find_full_range_haplotypes(
    node_info: &HashMap<String, NodeInfo>,
    all_sequences: &HashMap<usize, Vec<(Vec<String>, String, HashSet<String>)>>,
) -> HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)> {
    let mut best_paths = HashMap::new();
    for (hap_index, path_list) in all_sequences.iter() {
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

/// Assign every graph interval at least one node per haplotype.
///
/// This iterates over *all* intervals in the graph (not just the ones a haplotype's
/// reads already cover) so that no haplotype is left with a gap. Within each interval:
///   * an uncontested interval (a single node) is treated as homozygous backbone and
///     assigned to every haplotype, so it never blocks the haplotype-constrained DFS;
///   * a contested interval picks, per haplotype, the node with the highest read overlap,
///     falling back to the most-supported node when there is no overlap so that every
///     haplotype still gets a node there.
pub fn filter_haplotype_nodes(
    node_info: &HashMap<String, NodeInfo>,
    haplotype_nodes: &HashMap<usize, HashSet<String>>,
    haplotype_reads: &HashMap<usize, HashSet<String>>,
) -> HashMap<usize, HashSet<String>> {
    // All intervals across the whole graph, each with every parallel node it contains.
    let interval_all_nodes =
        find_parallele_nodes_from_nodelist(&node_info.keys().cloned().collect::<Vec<_>>());

    // The full set of haplotype indices we must cover in every interval.
    let mut hap_set: HashSet<usize> = HashSet::new();
    hap_set.extend(haplotype_reads.keys().cloned());
    hap_set.extend(haplotype_nodes.keys().cloned());

    let mut filtered_haplotype_nodes: HashMap<usize, HashSet<String>> = HashMap::new();
    for hap in hap_set.iter() {
        filtered_haplotype_nodes.entry(*hap).or_default();
    }

    let support_of = |node: &String| -> usize {
        node_info.get(node).map(|n| n.support_reads).unwrap_or(0)
    };

    for (interval, nodes_set) in interval_all_nodes.iter() {
        let nodes: Vec<String> = nodes_set.iter().cloned().collect();
        if nodes.is_empty() {
            continue;
        }

        // Uncontested interval: a single allele shared by every haplotype (homozygous backbone).
        if nodes.len() == 1 {
            let node = nodes[0].clone();
            for hap in hap_set.iter() {
                filtered_haplotype_nodes
                    .entry(*hap)
                    .or_default()
                    .insert(node.clone());
            }
            continue;
        }

        // Deterministic fallback: most-supported node (tie-break: lexicographically smaller id),
        // guaranteeing every haplotype is assigned a node even with no read overlap.
        let fallback_node = nodes
            .iter()
            .max_by(|a, b| support_of(a).cmp(&support_of(b)).then_with(|| b.cmp(a)))
            .cloned()
            .unwrap();

        for hap in hap_set.iter() {
            let read_list = haplotype_reads.get(hap).cloned().unwrap_or_default();

            let mut best_node: Option<String> = None;
            let mut best_intersection = 0usize;
            for node in nodes.iter() {
                let node_read_names = get_read_name_list(node_info, node.clone());
                let intersection_count = read_list.intersection(&node_read_names).count();
                let better = match &best_node {
                    None => intersection_count > 0,
                    Some(current) => {
                        if intersection_count > best_intersection {
                            true
                        } else if intersection_count == best_intersection {
                            // tie-break: higher support, then lexicographically smaller id
                            let ns = support_of(node);
                            let cs = support_of(current);
                            ns > cs || (ns == cs && node < current)
                        } else {
                            false
                        }
                    }
                };
                if better {
                    best_node = Some(node.clone());
                    best_intersection = intersection_count;
                }
            }

            let chosen = match best_node {
                Some(n) if best_intersection > 0 => n,
                _ => {
                    warn!(
                        "hap: {}, interval: {}, no read overlap; falling back to most-supported node {}",
                        hap, interval, fallback_node
                    );
                    fallback_node.clone()
                }
            };
            filtered_haplotype_nodes
                .entry(*hap)
                .or_default()
                .insert(chosen);
        }
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

        let haplotype_reads =
            assign_haplotype_reads(node_info, &heterozygous_nodes, hap_number);

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

    writeln!(file, "##fileformat=BED")?;
    writeln!(file, "##haplotype={}", haplotype_index + 1)?;
    writeln!(
        file,
        "#CHROM\tRef_start\tRef_end\tMod_rate\tAsm_start\tAsm_end\tMotif\tCoverage"
    )?;
    let mut methyl_info_vec = Vec::new();
    for ((ref_pos, asm_pos), (score, coverage)) in methyl_info.iter() {
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
    haplotype_number: usize,
    output_prefix: &PathBuf
) -> AnyhowResult<(HashMap<usize, (Vec<String>, String, HashSet<String>, usize, usize)>, HashMap<String, NodeInfo>, HashMap<String, Vec<String>>)> {
    let (node_info, edge_info) = load_graph(graph_filename).unwrap();

    let (_haplotype_reads, node_haplotype) =
        find_node_haplotype(&node_info, haplotype_number);
    let primary_haplotypes = find_best_haplotype_paths_dp(
        &node_info,
        &edge_info,
        &node_haplotype,
        haplotype_number,
    );
    info!("Best haplotype paths (DP): {}", primary_haplotypes.len());
    call_methylation(&node_info, primary_haplotypes.clone(), output_prefix);
    info!(
        "Haplotype specific methylation signals exported: {}",
        primary_haplotypes.len()
    );
    let output_filename = PathBuf::from(format!("{}.fasta", output_prefix.to_string_lossy()));
    let _ = write_graph_path_fasta(&primary_haplotypes, &output_filename);
    info!(
        "All sequences written to fasta: {}",
        output_prefix.to_str().unwrap()
    );
    Ok((primary_haplotypes, node_info.clone(), edge_info.clone()))
}

#[cfg(test)]
mod path_dp_tests {
    use super::{
        find_best_haplotype_paths_dp, find_full_range_haplotypes, enumerate_all_paths_with_haplotype,
        construct_sequences_from_haplotype_path, NodeInfo,
    };
    use std::collections::{HashMap, HashSet};

    fn test_node(id: &str, seq: &str, reads: &str, support: usize) -> NodeInfo {
        NodeInfo {
            seq: seq.to_string(),
            cigar: format!("{}M", seq.len()),
            support_reads: support,
            allele_frequency: "0.5".to_string(),
            read_names: reads.to_string(),
            methyl_info: HashMap::new(),
        }
    }

    #[test]
    fn dp_matches_longest_dfs_path_on_linear_graph() {
        let mut node_info = HashMap::new();
        node_info.insert(
            "H.chr1:1-10.0".to_string(),
            test_node("H.chr1:1-10.0", "AAAAAAAAAA", "r1,r2", 2),
        );
        node_info.insert(
            "H.chr1:10-20.0".to_string(),
            test_node("H.chr1:10-20.0", "CCCCCCCCCC", "r1,r2", 2),
        );
        node_info.insert(
            "H.chr1:20-30.0".to_string(),
            test_node("H.chr1:20-30.0", "GGGGGGGGGG", "r1,r2", 2),
        );

        let mut edge_info = HashMap::new();
        edge_info.insert(
            "H.chr1:1-10.0".to_string(),
            vec!["H.chr1:10-20.0".to_string()],
        );
        edge_info.insert(
            "H.chr1:10-20.0".to_string(),
            vec!["H.chr1:20-30.0".to_string()],
        );

        let mut node_haplotype = HashMap::new();
        for id in node_info.keys() {
            node_haplotype.insert(id.clone(), HashSet::from([0]));
        }

        let dp = find_best_haplotype_paths_dp(&node_info, &edge_info, &node_haplotype, 1);
        assert_eq!(dp.len(), 1);
        assert_eq!(dp[&0].0.len(), 3);
        assert_eq!(dp[&0].1.len(), 30);

        let all_paths = enumerate_all_paths_with_haplotype(
            &node_info,
            &edge_info,
            &node_haplotype,
            1,
        )
        .unwrap();
        let allseq = construct_sequences_from_haplotype_path(&node_info, &all_paths);
        let dfs_best = find_full_range_haplotypes(&node_info, &allseq);
        assert_eq!(dp[&0].1, dfs_best[&0].1);
    }
}