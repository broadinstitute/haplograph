use adjustp::{adjust, Procedure};
use anyhow::Result as AnyhowResult;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
use bio::io::fasta::Reader as FastaReader;
use bio::io::fastq;
use flate2::read::GzDecoder;
use indicatif::ProgressBar;
use log::{debug, warn};
use ndarray::s;
use ndarray::{Array1, Array2};
use rand::seq::SliceRandom;
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, IndexedReader, Read as BamRead};
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Normal;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use url::Url;
use std::path::PathBuf;
use std::io::Write;

pub fn reverse_complement(kmer: &str) -> String {
    kmer.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => panic!("Unexpected character: {}", c),
        })
        .collect()
}

pub fn write_fasta(
    all_sequences: &HashMap<String, String>,
    output_filename: &PathBuf,
) -> AnyhowResult<()> {
    if let Some(parent) = output_filename.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut file = File::create(output_filename)?;
    let chars_per_line = 60;
    for (header, sequence) in all_sequences.iter() {
       writeln!(file, ">{}", header)?;
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


pub fn gcs_gcloud_is_installed() -> bool {
    // Check if gcloud is installed on the PATH
    // Suppress stdout and stderr to prevent them from printing to the screen
    let mut cmd = std::process::Command::new("gcloud");
    cmd.arg("version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .is_ok()
}

pub fn gcs_authorize_data_access() {
    // Check if gcloud is installed on the PATH
    if !gcs_gcloud_is_installed() {
        panic!("gcloud is not installed on the PATH");
    }

    // Execute the command and capture the output
    let output = std::process::Command::new("gcloud")
        .args(["auth", "application-default", "print-access-token"])
        .output()
        .expect("Failed to execute command");

    if !output.status.success() {
        panic!("{}", String::from_utf8_lossy(&output.stderr));
    }

    // Decode the output and remove trailing newline
    let token = String::from_utf8(output.stdout)
        .expect("Failed to decode output")
        .trim_end()
        .to_string();

    // Set the environment variable
    unsafe {
        std::env::set_var("GCS_OAUTH_TOKEN", token);
    }
}

// Function to get a mapping between read group and sample name from a BAM header.
pub fn get_rg_to_sm_mapping(bam: &IndexedReader) -> HashMap<String, String> {
    let header = bam::Header::from_template(bam.header());

    let rg_sm_map: HashMap<String, String> = header
        .to_hashmap()
        .into_iter()
        .flat_map(|(_, records)| records)
        .filter(|record| record.contains_key("ID") && record.contains_key("SM"))
        .map(|record| (record["ID"].clone(), record["SM"].clone()))
        .collect();

    rg_sm_map
}

pub fn get_sm_name_from_rg(
    read: &bam::Record,
    rg_sm_map: &HashMap<String, String>,
) -> AnyhowResult<String> {
    let rg = read.aux(b"RG")?;

    if let Aux::String(v) = rg {
        if let Some(sm) = rg_sm_map.get(v) {
            Ok(sm.to_owned())
        } else {
            Err(anyhow::anyhow!(
                "Sample name not found for read group: {}",
                v
            ))
        }
    } else {
        Err(anyhow::anyhow!("Read group is not a string"))
    }
}

pub fn open_bam_file(alignment_bam: &String) -> IndexedReader {
    // Open BAM file
    if alignment_bam.starts_with("gs://") && std::env::var("GCS_OAUTH_TOKEN").is_err() {
        gcs_authorize_data_access();
    }
    let bam = if alignment_bam.starts_with("gs://") {
        let url = Url::parse(alignment_bam).unwrap();
        IndexedReader::from_url(&url).unwrap()
    } else {
        IndexedReader::from_path(alignment_bam).unwrap()
    };

    debug!("Successfully opened BAM file");

    bam
}

pub fn get_chromosome_ref_seq(reference_fa: &String, chromosome: &str) -> Vec<fastq::Record> {
    debug!("Opening FASTA file: {}", reference_fa);

    // Open FASTA file (supports both regular and gzipped files)
    let file = File::open(reference_fa).expect("Failed to open FASTA file");
    let reader: Box<dyn BufRead> = if reference_fa.ends_with(".gz") {
        let gz_decoder = GzDecoder::new(file);
        Box::new(BufReader::new(gz_decoder))
    } else {
        Box::new(BufReader::new(file))
    };
    let fasta_reader = FastaReader::new(reader);

    let mut reference_seqs = Vec::new();

    for result in fasta_reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        let sequence = String::from_utf8_lossy(record.seq()).to_string();

        // Check if this sequence matches the target chromosome
        if seq_id == chromosome {
            debug!(
                "Extracting region {} from chromosome {}",
                chromosome, seq_id
            );
            let chromosome_seq = sequence;
            let record_id = seq_id;
            let fastq_record = fastq::Record::with_attrs(
                &record_id,
                None,
                chromosome_seq.as_bytes(),
                vec![30; chromosome_seq.len()].as_slice(), // Default quality score
            );
            reference_seqs.push(fastq_record);
            debug!("Extracted reference sequence: {} bp", chromosome_seq.len());
        }
    }

    if reference_seqs.is_empty() {
        warn!(
            "No matching chromosome '{}' found in FASTA file",
            chromosome
        );
    }

    reference_seqs
}

pub fn get_all_ref_seq(reference_fa: &String) -> Vec<fastq::Record> {
    debug!("Opening FASTA file: {}", reference_fa);

    // Open FASTA file (supports both regular and gzipped files)
    let file = File::open(reference_fa).expect("Failed to open FASTA file");
    let reader: Box<dyn BufRead> = if reference_fa.ends_with(".gz") {
        let gz_decoder = GzDecoder::new(file);
        Box::new(BufReader::new(gz_decoder))
    } else {
        Box::new(BufReader::new(file))
    };
    let fasta_reader = FastaReader::new(reader);

    let mut reference_seqs = Vec::new();

    for result in fasta_reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        let sequence = String::from_utf8_lossy(record.seq()).to_string();

        let fastq_record = fastq::Record::with_attrs(
            &seq_id,
            None,
            sequence.as_bytes(),
            vec![30; sequence.len()].as_slice(), // Default quality score
        );
        reference_seqs.push(fastq_record);
        // info!("Extracted {} reference sequences", reference_seqs.len());
    }

    reference_seqs
}

pub fn get_ref_seq_from_chromosome(
    reference_fa: &String,
    chromosome: &str,
) -> (Vec<fastq::Record>, Vec<fastq::Record>) {
    let reference_seqs = get_all_ref_seq(reference_fa);
    let reference_seqs_chromosome = reference_seqs
        .iter()
        .filter(|r| r.id() == chromosome).cloned()
        .collect::<Vec<fastq::Record>>();

    (reference_seqs, reference_seqs_chromosome)
}

pub fn split_locus(locus: String) -> (String, usize, usize) {
    let parts: Vec<&str> = locus.split(':').collect();
    let chromosome = parts[0].to_string();
    let start: usize = parts[1]
        .split("-")
        .collect::<Vec<&str>>()
        .first()
        .unwrap()
        .parse()
        .unwrap();
    let end: usize = parts[1]
        .split("-")
        .collect::<Vec<&str>>()
        .last()
        .unwrap()
        .parse()
        .unwrap();
    (chromosome, start, end)
}

pub fn mask_ns(seq: &str) -> String {
    seq.chars()
        .map(|c| {
            let upper_c = c.to_ascii_uppercase();
            if !['A', 'G', 'C', 'T'].contains(&upper_c) {
                c.to_ascii_lowercase()
            } else {
                upper_c
            }
        })
        .collect()
}

pub fn alignment_to_cigar(operations: &[AlignmentOperation]) -> String {
    let mut cigar: Vec<(usize, char)> = Vec::new();

    for op in operations {
        let cigar_op = match op {
            AlignmentOperation::Match => '=',
            AlignmentOperation::Subst => 'X',
            AlignmentOperation::Del => 'D',
            AlignmentOperation::Ins => 'I',
            AlignmentOperation::Xclip(_) => 'S',
            AlignmentOperation::Yclip(_) => 'S',
        };

        if !cigar.is_empty() && cigar.last().unwrap().1 == cigar_op {
            cigar.last_mut().unwrap().0 += 1;
        } else {
            cigar.push((1, cigar_op));
        }
    }

    cigar
        .iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect()
}

pub fn gap_open_aligner(reference: &str, sequence: &str) -> String {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // Create an aligner with the same scoring parameters
    let mut aligner = Aligner::with_capacity(sequence.len(), reference.len(), -5, -1, &score); // match_score=0, mismatch_score=-6, gap_open=-5, gap_extend=-3

    // Perform the alignment
    let alignment = aligner.global(sequence.as_bytes(), reference.as_bytes());

    // Get the aligned sequences
    
    // println!("{:?}", cigar);

    alignment_to_cigar(&alignment.operations)
}

/// Find overlapping reads between two read vectors
pub fn find_overlapping_reads(read_vector1: &[String], read_vector2: &[String]) -> Vec<String> {
    let set1: std::collections::HashSet<String> = read_vector1
        .iter()
        .map(|s| s.split('|').next().unwrap().to_string())
        .collect();
    let set2: std::collections::HashSet<String> = read_vector2
        .iter()
        .map(|s| s.split('|').next().unwrap().to_string())
        .collect();

    set1.intersection(&set2).cloned()
        .collect()
}

pub fn import_bed(bed_file: &String) -> Vec<(String, usize, usize)> {
    let mut bed_list = Vec::new();
    if bed_file.ends_with(".gz") {
        let file = File::open(bed_file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            bed_list.push((
                line.split('\t').next().unwrap().to_string(),
                line.split('\t').nth(1).unwrap().parse().unwrap(),
                line.split('\t').nth(2).unwrap().parse().unwrap(),
            ));
        }
    } else {
        let file = File::open(bed_file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            bed_list.push((
                line.split('\t').next().unwrap().to_string(),
                line.split('\t').nth(1).unwrap().parse().unwrap(),
                line.split('\t').nth(2).unwrap().parse().unwrap(),
            ));
        }
    }
    bed_list
}

pub fn process_cigar(cigar: &str) -> String {
    let mut out = String::new();
    let mut n = 0;

    for symbol in cigar.chars() {
        if symbol.is_ascii_digit() {
            n = 10 * n + symbol.to_digit(10).unwrap() as usize;
        } else {
            if n == 0 {
                out.push(symbol);
            } else {
                out.push_str(&symbol.to_string().repeat(n));
            }
            n = 0;
        }
    }

    out
}

pub fn combine_cigar(cigar: &str) -> String {
    if cigar.is_empty() {
        return String::new();
    }

    // Convert to Vec<char> for efficient indexing
    let mut chars: Vec<char> = cigar.chars().collect();
    chars.push('$'); // Add sentinel character to handle the last group

    let mut out = String::new();
    let mut start = 0;

    for i in 1..chars.len() {
        if chars[i - 1] != chars[i] {
            let length = i - start;
            out.push_str(&format!("{}{}", length, chars[i - 1]));
            start = i;
        }
    }

    out
}

fn jaccard_distance(vector1: &[bool], vector2: &[bool]) -> f64 {
    assert_eq!(
        vector1.len(),
        vector2.len(),
        "Vectors must have the same length"
    );

    let mut intersection_count = 0;
    let mut union_count = 0;

    for (a, b) in vector1.iter().zip(vector2.iter()) {
        if *a && *b {
            intersection_count += 1;
        }
        if *a || *b {
            union_count += 1;
        }
    }

    if union_count == 0 {
        return 0.0; // Both vectors are all zeros
    }

    1.0 - (intersection_count as f64 / union_count as f64)
}

/// Generate a null distribution through permutation testing
fn get_null_distribution(
    records: &Vec<String>,
    matrix: &Array2<f64>,
    permutation_round: usize,
) -> Vec<f64> {
    let bar = ProgressBar::new(permutation_round as u64);
    let summary_statistics = (0..permutation_round)
        .into_par_iter()
        .flat_map(|_| {
            bar.inc(1);
            let mut local_stats = Vec::new();
            let mut rng = rand::rng();

            for (i, index) in records.iter().enumerate() {
                let vector = matrix.slice(s![i, ..]);

                // Create a shuffled copy of the vector
                let vector_data: Vec<f64> = vector.iter().copied().collect();
                let mut shuffled_data = vector_data.clone();
                shuffled_data.shuffle(&mut rng);
                let shuffled = Array1::from(shuffled_data);

                let mut all_coefficients = Vec::new();

                for (j, other_index) in records.iter().enumerate() {
                    if index == other_index {
                        continue;
                    }

                    let other_vector = matrix.slice(s![j, ..]);

                    let binary_vector: Vec<bool> = shuffled.iter().map(|&x| x > 0.5).collect();
                    let binary_other: Vec<bool> = other_vector.iter().map(|&x| x > 0.5).collect();

                    // Calculate Jaccard distance
                    let coor = (1.0 - jaccard_distance(&binary_vector, &binary_other)).abs();
                    all_coefficients.push(coor);
                }

                local_stats.push(all_coefficients.iter().sum());
            }
            local_stats
        })
        .collect::<Vec<f64>>();
    bar.finish();
    summary_statistics
}

/// Calculate statistics for observed data
fn calculate_observation_statistics(
    recordlist: &Vec<String>,
    index: usize,
    matrix: &Array2<f64>,
) -> f64 {
    let vector = &matrix.slice(s![index, ..]);
    let mut all_coefficients = Vec::new();

    for (i, other_index) in recordlist.iter().enumerate() {
        if i == index {
            continue;
        }

        let other_vector = &matrix.slice(s![i, ..]);

        // Convert arrays to binary vectors before calculating Jaccard distance
        let binary_vector: Vec<bool> = vector.iter().map(|&x| x > 0.5).collect();
        let binary_other: Vec<bool> = other_vector.iter().map(|&x| x > 0.5).collect();

        // Calculate Jaccard distance
        let coor = (1.0 - jaccard_distance(&binary_vector, &binary_other)).abs();
        all_coefficients.push(coor);
    }

    all_coefficients.iter().sum()
}

/// Calculate p-value using z-score approach
fn calculate_p_value(statistics: &[f64], observation: f64) -> f64 {
    let n = statistics.len() as f64;

    // Calculate mean
    let mu = statistics.iter().sum::<f64>() / n;

    // Calculate standard deviation
    let variance = statistics.iter().map(|&x| (x - mu).powi(2)).sum::<f64>() / n;
    let sigma = variance.sqrt();
    // println!("{:?}, {}", statistics, observation);

    let z_score = (observation - mu) / sigma;

    // Calculate p-value using normal distribution CDF
    let normal = Normal::new(0.0, 1.0).unwrap();
    1.0 - normal.cdf(z_score)
}

pub fn permutation_test(
    matrix: &Array2<f64>,
    p_value_threshold: f64,
    permutation_round: usize,
    node_list: Vec<String>,
) -> Vec<String> {
    let bar = ProgressBar::new(node_list.len() as u64);
    let statistics = get_null_distribution(&node_list, matrix, permutation_round);
    // Replace par_iter().enumerate() with this pattern
    let (indices, collected_values): (Vec<_>, Vec<_>) = (0..node_list.len())
        .into_par_iter()
        .map(|i| {
            bar.inc(1);
            let nodename = &node_list[i];
            let observation = calculate_observation_statistics(&node_list, i, matrix);
            let p_value = calculate_p_value(&statistics, observation);
            ((nodename.clone()), Some((p_value, nodename.clone())))
        })
        .unzip();

    let mut raw_p_values = Vec::new();
    let mut test_index = Vec::new();

    for item in collected_values.into_iter().flatten() {
        let (p_value, node_name) = item;
        if !p_value.is_nan() {
            raw_p_values.push(p_value);
            test_index.push(node_name);
        }
    }

    // adjust pvalues, create excluded_index list
    let mut excluded_index = Vec::new();
    // println!("raw_p_values: {:?}", raw_p_values);
    if !raw_p_values.is_empty() {
        let qvalues = adjust(&raw_p_values, Procedure::BenjaminiHochberg);
        for (qi, q_value) in qvalues.iter().enumerate() {
            let test_index_value = &test_index[qi];
            if q_value > &p_value_threshold {
                excluded_index.push(test_index_value);
                debug!(
                    "excluded_index: {:?}, q_value: {:?}",
                    test_index_value, q_value
                );
            }
        }
    }

    bar.finish();

    // filter variants
    let mut index_list = Vec::new();
    let mut f_node: Vec<String> = Vec::new();
    // get index list and var_list
    for (r, rindex) in node_list.iter().enumerate() {
        if !excluded_index.contains(&rindex) {
            index_list.push(r);
            f_node.push(rindex.clone());
        }
    }

    f_node
}
