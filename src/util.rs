use std::collections::{HashMap, HashSet};
use url::Url;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use log::{info, warn};
use anyhow::{Result as AnyhowResult};
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
use bio::io::fasta::Reader as FastaReader;
use bio::io::fastq;

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
    std::env::set_var("GCS_OAUTH_TOKEN", token);
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

pub fn get_sm_name_from_rg(read: &bam::Record, rg_sm_map: &HashMap<String, String>) -> AnyhowResult<String> {
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
    let mut bam = if alignment_bam.starts_with("gs://") {
        let url = Url::parse(&alignment_bam).unwrap();
        IndexedReader::from_url(&url).unwrap()
    } else {
        IndexedReader::from_path(&alignment_bam).unwrap()
    };
    
    info!("Successfully opened BAM file");

    bam
}

pub fn get_chromosome_ref_seq(reference_fa: &String, chromosome: &str) -> Vec<fastq::Record> {
    
    
    info!("Opening FASTA file: {}", reference_fa);
    
    // Open FASTA file
    let mut fasta_reader = FastaReader::from_file(&reference_fa)
        .expect("Failed to open FASTA file");

    let mut reference_seqs = Vec::new();

    for result in fasta_reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        let sequence = String::from_utf8_lossy(record.seq()).to_string();

        // Check if this sequence matches the target chromosome
        if seq_id == chromosome {
            info!("Extracting region {} from chromosome {}", 
                chromosome, seq_id);
            let chromosome_seq = sequence;
            let record_id = seq_id;
            let fastq_record = fastq::Record::with_attrs(
                &record_id,
                None,
                chromosome_seq.as_bytes(),
                vec![30; chromosome_seq.len()].as_slice() // Default quality score
            );
            reference_seqs.push(fastq_record);
            info!("Extracted reference sequence: {} bp", chromosome_seq.len());
            
        }
    }
    
    if reference_seqs.is_empty() {
        warn!("No matching chromosome '{}' found in FASTA file", chromosome);
    }
    
    reference_seqs
}

pub fn get_all_ref_seq(reference_fa: &String) -> Vec<fastq::Record> {
    
    info!("Opening FASTA file: {}", reference_fa);
    
    // Open FASTA file
    let mut fasta_reader = FastaReader::from_file(&reference_fa)
        .expect("Failed to open FASTA file");

    let mut reference_seqs = Vec::new();

    for result in fasta_reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        let sequence = String::from_utf8_lossy(record.seq()).to_string();

        let fastq_record = fastq::Record::with_attrs(
            &seq_id,
            None,
            sequence.as_bytes(),
            vec![30; sequence.len()].as_slice() // Default quality score
        );
        reference_seqs.push(fastq_record);
        info!("Extracted {} reference sequences", reference_seqs.len());
            
        
    }
    
    reference_seqs
}

pub fn split_locus(locus: String) -> (String, usize, usize) {
    let parts: Vec<&str> = locus.split(':').collect();
    let chromosome = parts[0].to_string();
    let start: usize = parts[1].split("-").collect::<Vec<&str>>().first().unwrap().parse().unwrap();
    let end: usize = parts[1].split("-").collect::<Vec<&str>>().last().unwrap().parse().unwrap();
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
    let cigar = alignment_to_cigar(&alignment.operations);
    // println!("{:?}", cigar);

    cigar
}

/// Find overlapping reads between two read vectors
pub fn find_overlapping_reads(read_vector1: &[String], read_vector2: &[String]) -> Vec<String> {
    let set1: std::collections::HashSet<String> = read_vector1.iter().map(|s| s.split('|').next().unwrap().to_string()).collect();
    let set2: std::collections::HashSet<String> = read_vector2.iter().map(|s| s.split('|').next().unwrap().to_string()).collect();
    
    set1.intersection(&set2)
        .map(|read_name| read_name.clone())
        .collect()
}

