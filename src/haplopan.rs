use crate::util;
use anyhow::Result;
use bio::io::fasta::Reader as FastaReader;
use flate2::read::GzDecoder;
use log::info;
use ndarray::Array2;
use rayon::prelude::*;
use rust_htslib::bam::Read as BamRead;
use std::collections::{HashMap, HashSet};
use std::f64;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

fn encode_base(c: u8) -> Option<u64> {
    match c.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

// Rolling canonical k-mer counts using 2-bit encoding (k ≤ 32).
// canonical = min(forward, reverse_complement) to collapse strand ambiguity.
// Non-ACGT bases reset the rolling window.
pub fn count_kmer_rolling(sequence: &[u8], rollingkmer_list: &[usize]) -> HashMap<u64, i32> {
    let mut counts: HashMap<u64, i32> = HashMap::new();
    for &k in rollingkmer_list {
        assert!(k <= 32, "k={} exceeds maximum supported k of 32", k);
        let mask: u64 = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        let rc_shift = 2 * (k - 1);
        let mut fwd: u64 = 0;
        let mut rev: u64 = 0;
        let mut valid: usize = 0;
        for &base in sequence {
            if let Some(b) = encode_base(base) {
                fwd = ((fwd << 2) | b) & mask;
                rev = (rev >> 2) | ((b ^ 3) << rc_shift);
                valid += 1;
            } else {
                fwd = 0;
                rev = 0;
                valid = 0;
            }
            if valid >= k {
                *counts.entry(fwd.min(rev)).or_insert(0) += 1;
            }
        }
    }
    counts
}

#[allow(clippy::type_complexity)]
pub fn load_pangenome(
    pangenome_path: &PathBuf,
    rollingkmer_list: &[usize],
) -> Result<(HashMap<String, HashMap<u64, i32>>, HashMap<String, String>)> {
    info!("Open Pangenome File: {}", pangenome_path.display());
    let file = File::open(pangenome_path).expect("Failed to open FASTA file");
    let reader: Box<dyn BufRead> = if pangenome_path.ends_with(".gz") {
        let gz_decoder = GzDecoder::new(file);
        Box::new(BufReader::new(gz_decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    let fasta_reader = FastaReader::new(reader);
    let records: Vec<_> = fasta_reader
        .records()
        .map(|result| result.expect("Failed to read FASTA record"))
        .collect();

    let results: Vec<_> = records
        .par_iter()
        .map(|record| {
            let seq_id = record.id().to_string();
            let seq_description = record.desc().unwrap_or("").to_string();
            let allele_id = seq_description.split(' ').next().unwrap_or("").to_string();
            let kmer_count = count_kmer_rolling(record.seq(), rollingkmer_list);
            (seq_id, allele_id, kmer_count)
        })
        .collect();

    let mut kmer_dict = HashMap::new();
    let mut seq_info = HashMap::new();
    for (seq_id, allele_id, kmer_count) in results {
        seq_info.insert(seq_id.clone(), allele_id);
        kmer_dict.insert(seq_id, kmer_count);
    }

    info!("Pangenome Extraction Completed");
    Ok((kmer_dict, seq_info))
}

#[allow(clippy::type_complexity)]
pub fn count_kmer_from_read(
    pangenome_kmer_dict: &HashMap<String, HashMap<u64, i32>>,
    read_path: &Path,
    locus: &str,
    rollingkmer_list: &[usize],
) -> Result<(HashMap<u64, i32>, HashMap<String, (String, Vec<u8>)>)> {
    let mut read_kmer_dict: HashMap<u64, i32> = HashMap::new();
    let mut read_seq_dict: HashMap<String, (String, Vec<u8>)> = HashMap::new();
    let pangenome_kmer_set: HashSet<u64> = pangenome_kmer_dict
        .values()
        .flat_map(|kd| kd.keys().copied())
        .collect();
    info!("Pangenome K-mer Set: {:?}", pangenome_kmer_set.len());
    // info!(
    //     "Read File: {}, {:?}",
    //     read_path.display(),
    //     read_path.display().to_string().ends_with(".bam")
    // );
    let (chromosome, start, end) = util::split_locus(locus.to_string());

    if read_path.display().to_string().ends_with(".bam") {
        let mut bam = util::open_bam_file(&read_path.display().to_string());
        let mut read_set = HashSet::new();
        let header = bam.header().to_owned();
        let tid = header
            .tid(chromosome.as_bytes())
            .ok_or_else(|| anyhow::anyhow!("Chromosome {} not found in BAM header", chromosome))?;

        bam.fetch((tid, start as i64, end as i64))
            .expect("Failed to fetch BAM records");
        for read in bam.records() {
            let record = read.expect("Failed to read BAM record");
            let read_id = String::from_utf8_lossy(record.qname()).into_owned();
            if read_set.contains(&read_id) {
                continue;
            }
            read_set.insert(read_id.clone());
            let seq_bytes = record.seq().as_bytes();
            let read_seq = String::from_utf8_lossy(&seq_bytes).into_owned();
            let read_qual = record.qual();
            let local_kmer_dict = count_kmer_rolling(&seq_bytes, rollingkmer_list);
            let mut total_count = 0;
            for (&kmer, &count) in local_kmer_dict.iter() {
                if !pangenome_kmer_set.contains(&kmer) {
                    continue;
                }
                total_count += count;
                *read_kmer_dict.entry(kmer).or_insert(0) += count;
            }
            if total_count > 0 {
                read_seq_dict.insert(read_id, (read_seq, read_qual.to_vec()));
            }
        }
    } else if read_path.display().to_string().ends_with(".fastq")
        || read_path.display().to_string().ends_with(".fq")
        || read_path.display().to_string().ends_with(".fa")
        || read_path.display().to_string().ends_with(".fasta")
    {
        return Err(anyhow::anyhow!(
            "file format: {} UnderConstruction",
            read_path.display()
        ));
    } else {
        return Err(anyhow::anyhow!(
            "Unsupported file format: {}",
            read_path.display()
        ));
    };

    Ok((read_kmer_dict, read_seq_dict))
}

fn calculate_matrix(
    read_kmer_dict: &HashMap<u64, i32>,
    pangenome_kmer_dict: &HashMap<String, HashMap<u64, i32>>,
    sample_name: &String,
) -> Result<(Array2<f64>, Vec<u64>, Vec<String>)> {
    let read_kmer_set: HashSet<u64> = read_kmer_dict.keys().copied().collect();
    info!("Read K-mer Set: {:?}", read_kmer_set.len());
    let mut mrow_index: Vec<u64> = read_kmer_set.into_iter().collect();
    mrow_index.sort_unstable();
    let mrow_index_dict: HashMap<u64, usize> = mrow_index
        .iter()
        .enumerate()
        .map(|(i, &kmer)| (kmer, i))
        .collect();

    let mut reference_index: Vec<String> = pangenome_kmer_dict.keys().cloned().collect();
    reference_index.sort();
    reference_index.push(sample_name.clone());
    let mcol_index_dict: HashMap<String, usize> = reference_index
        .iter()
        .enumerate()
        .map(|(i, hap)| (hap.clone(), i))
        .collect();
    info!(
        "Matrix Shape: {:?}, {:?}",
        reference_index.len(),
        mrow_index.len()
    );

    let mut matrix = Array2::<f64>::zeros((mrow_index.len(), reference_index.len()));
    for haplotype in reference_index.iter() {
        let kmer_dict: &HashMap<u64, i32> = if haplotype == sample_name {
            read_kmer_dict
        } else {
            pangenome_kmer_dict.get(haplotype).unwrap()
        };
        for (&kmer, &count) in kmer_dict.iter() {
            if let Some(&row_index) = mrow_index_dict.get(&kmer) {
                let col_index = mcol_index_dict[haplotype];
                matrix[[row_index, col_index]] = count as f64;
            }
        }
    }

    Ok((matrix, mrow_index, reference_index))
}

fn select_haplotypes(
    x: &Array2<f64>,
    entities: &[String],
    sample_id: &String,
    ploidy: usize,
) -> (Vec<String>, f64) {
    let n = entities.len();
    let (nrows, _) = x.dim();
    let t_idx = entities
        .iter()
        .position(|s| s == sample_id)
        .expect("sample not found in entities");

    let hap_indices: Vec<usize> = (0..n).filter(|&i| i != t_idx).collect();

    let sample_total: f64 = (0..nrows).map(|r| x[[r, t_idx]]).sum();
    let mean_hap_total: f64 = {
        let sum: f64 = hap_indices
            .iter()
            .map(|&j| (0..nrows).map(|r| x[[r, j]]).sum::<f64>())
            .sum();
        if hap_indices.is_empty() {
            1.0
        } else {
            sum / hap_indices.len() as f64
        }
    };
    let d = if mean_hap_total == 0.0 || ploidy == 0 {
        1.0
    } else {
        sample_total / (ploidy as f64 * mean_hap_total)
    };

    let s_norm: Vec<f64> = (0..nrows).map(|r| x[[r, t_idx]] / d).collect();
    let s_norm_sq: f64 = s_norm.iter().map(|&v| v * v).sum();

    let mut remaining = s_norm.clone();
    let mut selected: Vec<String> = Vec::with_capacity(ploidy);

    for _ in 0..ploidy {
        let scores: Vec<f64> = hap_indices
            .iter()
            .copied()
            .map(|j| {
                (0..nrows)
                    .map(|r| {
                        let diff = remaining[r] - x[[r, j]];
                        diff * diff
                    })
                    .sum::<f64>()
            })
            .collect();
        let best_local = scores
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| hap_indices[i]);

        if let Some(best) = best_local {
            for r in 0..nrows {
                remaining[r] -= x[[r, best]];
            }
            selected.push(entities[best].clone());
        }
    }

    let residual_sq: f64 = remaining.iter().map(|&v| v * v).sum();
    let distance = if s_norm_sq == 0.0 {
        0.0
    } else {
        residual_sq / s_norm_sq
    };

    (selected, distance)
}

fn run_cmd(mut cmd: Command) -> anyhow::Result<()> {
    let prog = cmd.get_program().to_string_lossy().into_owned();
    let out = cmd.output()?;
    if !out.status.success() {
        return Err(anyhow::anyhow!(
            "{} failed (exit {}): {}",
            prog,
            out.status,
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }
    Ok(())
}

// Align `fastq_path` reads against `ref_fa`, producing `{stem}.sorted.bam` + `.bai`.
fn align_and_sort(ref_fa: &Path, fastq_path: &Path, out_stem: &Path) -> anyhow::Result<()> {
    let sam_path = out_stem.with_extension("sam");
    let unsorted_bam = out_stem.with_extension("bam");
    let sorted_bam = PathBuf::from(format!("{}.sorted.bam", out_stem.display()));

    let mut minimap2 = Command::new("minimap2");
    minimap2
        .args(["-ax", "asm5", "--cs", "--end-bonus=10"])
        .arg(ref_fa)
        .arg(fastq_path)
        .args(["-o", sam_path.to_str().unwrap()]);
    run_cmd(minimap2)?;

    let mut view = Command::new("samtools");
    view.args(["view", "-bS"])
        .arg(&sam_path)
        .args(["-o"])
        .arg(&unsorted_bam);
    run_cmd(view)?;

    let mut sort = Command::new("samtools");
    sort.arg("sort")
        .arg(&unsorted_bam)
        .args(["-o"])
        .arg(&sorted_bam);
    run_cmd(sort)?;

    let mut index = Command::new("samtools");
    index.arg("index").arg(&sorted_bam);
    run_cmd(index)?;

    std::fs::remove_file(&sam_path).ok();
    std::fs::remove_file(&unsorted_bam).ok();
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn start(
    read_path: &Path,
    pangenome_path: &PathBuf,
    locus: &str,
    rollingkmer_list: &[usize],
    output_prefix: &String,
    _data_technology: &str,
    sample_id: &String,
    ploidy: usize,
) -> Result<(Vec<String>, f64)> {
    info!("Loading Pangenome File: {}", pangenome_path.display());
    let (pangenome_kmer_dict, _seq_info) = load_pangenome(pangenome_path, rollingkmer_list)?;
    info!("Counting K-mers from Read File: {}", read_path.display());
    let (read_kmer_dict, read_seq_dict) =
        count_kmer_from_read(&pangenome_kmer_dict, read_path, locus, rollingkmer_list)?;
    info!("K-mer Counting Completed");
    let (matrix, _mrow_index, reference_index) =
        calculate_matrix(&read_kmer_dict, &pangenome_kmer_dict, sample_id)?;

    info!("Calculating Best Haplotypes");
    let (selected, distance) = select_haplotypes(&matrix, &reference_index, sample_id, ploidy);
    info!(
        "Selected haplotypes: {:?}, distance: {:.4}",
        selected, distance
    );

    let mut seen = std::collections::HashSet::new();
    let unique_selected: Vec<String> = selected
        .iter()
        .filter(|h| seen.insert(h.as_str()))
        .cloned()
        .collect();

    let mut final_sequences_dict = HashMap::new();
    let fasta_reader = FastaReader::from_file(pangenome_path)?;
    for record in fasta_reader.records() {
        let record = record.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        if unique_selected.contains(&seq_id) {
            final_sequences_dict.insert(
                seq_id.replace(":", "_"),
                String::from_utf8_lossy(record.seq()).to_string(),
            );
        }
    }
    info!("Writing Best Haplotypes to File");
    let tmp_fasta_path = PathBuf::from(format!("{}.tmp.ref.fasta", output_prefix));
    util::write_fasta(&final_sequences_dict, &tmp_fasta_path)?;

    // Write locus reads already in memory to FASTQ — avoids re-reading the whole BAM.
    let tmp_fastq_path = PathBuf::from(format!("{}.tmp.fastq", output_prefix));
    {
        use std::io::Write;
        let mut fq = std::io::BufWriter::new(std::fs::File::create(&tmp_fastq_path)?);
        for (read_id, (seq, qual)) in read_seq_dict.iter() {
            writeln!(fq, "@{}", read_id)?;
            writeln!(fq, "{}", seq)?;
            writeln!(fq, "+")?;
            // rust-htslib qual() returns raw Phred (0-based); FASTQ needs Phred+33
            let qual_ascii: Vec<u8> = qual.iter().map(|&q| q + 33).collect();
            fq.write_all(&qual_ascii)?;
            writeln!(fq)?;
        }
    }
    info!(
        "Realigning {} locus reads to selected haplotypes",
        read_seq_dict.len()
    );
    let out_stem = PathBuf::from(format!("{}.tmp", output_prefix));
    align_and_sort(&tmp_fasta_path, &tmp_fastq_path, &out_stem)?;
    std::fs::remove_file(&tmp_fastq_path).ok();

    Ok((selected, distance))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    // ── select_haplotypes tests ──────────────────────────────────────────────

    #[test]
    fn test_select_haplotypes_heterozygous() {
        // hap0=[2,0,1,0], hap1=[0,3,0,1], sample=[2,3,1,1] = hap0+hap1
        // D = 7 / (2 * 3.5) = 1.0  =>  s_norm = sample
        let data = vec![2.0, 0.0, 2.0, 0.0, 3.0, 3.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0];
        let x = Array2::from_shape_vec((4, 3), data).unwrap();
        let entities = vec!["hap0".to_string(), "hap1".to_string(), "sample".to_string()];
        let (selected, distance) = select_haplotypes(&x, &entities, &"sample".to_string(), 2);
        assert_eq!(selected.len(), 2);
        assert!(
            selected.contains(&"hap0".to_string()),
            "hap0 not in {:?}",
            selected
        );
        assert!(
            selected.contains(&"hap1".to_string()),
            "hap1 not in {:?}",
            selected
        );
        assert!(
            distance < 0.01,
            "expected near-zero distance, got {:.4}",
            distance
        );
    }

    #[test]
    fn test_select_haplotypes_homozygous() {
        // hap0=[2,1,3], hap1=[0,5,0] (very different), sample=[4,2,6] = 2*hap0
        let data = vec![2.0, 0.0, 4.0, 1.0, 5.0, 2.0, 3.0, 0.0, 6.0];
        let x = Array2::from_shape_vec((3, 3), data).unwrap();
        let entities = vec!["hap0".to_string(), "hap1".to_string(), "sample".to_string()];
        let (selected, _distance) = select_haplotypes(&x, &entities, &"sample".to_string(), 2);
        assert_eq!(selected.len(), 2);
        assert_eq!(selected[0], "hap0");
        assert_eq!(selected[1], "hap0");
    }

    #[test]
    fn test_select_haplotypes_haploid() {
        // hap0=[2,0], hap1=[0,3], sample=[2,0] matches hap0 only
        let data = vec![2.0, 0.0, 2.0, 0.0, 3.0, 0.0];
        let x = Array2::from_shape_vec((2, 3), data).unwrap();
        let entities = vec!["hap0".to_string(), "hap1".to_string(), "sample".to_string()];
        let (selected, _distance) = select_haplotypes(&x, &entities, &"sample".to_string(), 1);
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0], "hap0");
    }

    // ── count_kmer_rolling tests ─────────────────────────────────────────────

    #[test]
    fn test_count_kmer_rolling_canonical() {
        // "ACG" and its reverse complement "CGT" must map to the same canonical u64 key.
        // fwd("ACG")=6, rev("ACG")=27 → canonical=6
        // fwd("CGT")=27, rev("CGT")=6 → canonical=6
        let counts_acg = count_kmer_rolling(b"ACG", &[3]);
        let counts_cgt = count_kmer_rolling(b"CGT", &[3]);
        assert_eq!(counts_acg.len(), 1);
        assert_eq!(counts_cgt.len(), 1);
        let key_acg = *counts_acg.keys().next().unwrap();
        let key_cgt = *counts_cgt.keys().next().unwrap();
        assert_eq!(key_acg, key_cgt, "ACG and RC CGT must share canonical key");
    }

    #[test]
    fn test_count_kmer_rolling_counts() {
        // "AAAA" with k=2: three overlapping AA 2-mers.
        // AA encodes fwd=0, RC(AA)=TT=0b1111=15 → canonical=min(0,15)=0.
        let counts = count_kmer_rolling(b"AAAA", &[2]);
        assert_eq!(counts.len(), 1);
        let (&key, &count) = counts.iter().next().unwrap();
        assert_eq!(key, 0u64, "canonical of AA should be 0");
        assert_eq!(count, 3, "AAAA has 3 overlapping 2-mers");
    }

    #[test]
    fn test_count_kmer_rolling_n_resets_window() {
        // "ACNGT" with k=2: the N resets the window, so CG is never emitted.
        // "ACGT"  with k=2: AC, CG, GT are all emitted (3 total k-mer instances).
        let counts_n = count_kmer_rolling(b"ACNGT", &[2]);
        let counts_full = count_kmer_rolling(b"ACGT", &[2]);
        let total_n: i32 = counts_n.values().sum();
        let total_full: i32 = counts_full.values().sum();
        assert!(total_n < total_full, "N should reduce k-mer count");
    }
}
