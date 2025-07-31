use std::collections::HashMap;
use std::io::Result;
use std::path::Path;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::io::{Read, Write};
use std::io::{Seek, SeekFrom};
use bio::io::fastq;
use rust_htslib::bam::{self, FetchDefinition, IndexedReader, Read as BamRead, record::Aux};
use rust_htslib::faidx::Reader;
use anyhow::{Result as AnyhowResult, Context};
use log::info;
use indicatif::{ProgressBar, ProgressStyle};


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
fn get_rg_to_sm_mapping(bam: &IndexedReader) -> HashMap<String, String> {
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

fn get_sm_name_from_rg(read: &bam::Record, rg_sm_map: &HashMap<String, String>) -> AnyhowResult<String> {
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


/// Extract reads from a single locus with additional filtering options
/// 
/// # Arguments
/// * `bam` - Mutable reference to indexed BAM reader
/// * `chr` - Chromosome name
/// * `start` - Start position (0-based)
/// * `end` - End position (0-based)
/// * `min_mapping_quality` - Minimum mapping quality filter (optional)
/// 
/// # Returns
/// * `Result<Vec<fastq::Record>>` - Vector of extracted reads
/// 
pub fn extract_reads_at_single_locus(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    min_mapping_quality: Option<u8>,
) -> AnyhowResult<Vec<fastq::Record>> {
    let rg_sm_map = get_rg_to_sm_mapping(bam);
    let mut bmap = HashMap::new();
    let mut read_spans: HashMap<String, (u64, u64)> = HashMap::new(); // Track read alignment spans
    
    let _ = bam.fetch((chr.as_bytes(), start, end));
    
    // Create progress bar
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap()
    );
    pb.set_message("Processing pileup...");
    
    for p in bam.pileup() {
        let pileup = p?;
        pb.set_message(format!("Processing position {}", pileup.pos()));

        if start <= (pileup.pos() as u64) && (pileup.pos() as u64) < end {
            for (i, alignment) in pileup.alignments().enumerate() {
                let record = alignment.record();
                
                // Apply mapping quality filter
                if let Some(min_qual) = min_mapping_quality {
                    if record.mapq() < min_qual {
                        continue;
                    }
                }
                
                let qname = String::from_utf8_lossy(record.qname()).into_owned();
                let sm = match get_sm_name_from_rg(&record, &rg_sm_map) {
                    Ok(a) => a,
                    Err(_) => String::from("unknown"),
                };

                let is_secondary = record.is_secondary();
                let is_supplementary = record.is_supplementary();
                let locus_key = format!("{}:{}-{}", chr, start, end);
                let seq_name = format!("{qname}|{locus_key}|{sm}");

                // Track alignment span for this read
                let reference_start = record.pos() as u64;
                assert!(record.cigar().pos() == reference_start as i64);
                let reference_end = record.cigar().end_pos() as u64; 
                let read_span = read_spans.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
                read_span.0 = read_span.0.min(reference_start);
                read_span.1 = read_span.1.max(reference_end);

                if !bmap.contains_key(&seq_name) {
                    bmap.insert(seq_name.clone(), (String::new(), Vec::new()));
                }

                if !alignment.is_del() && !alignment.is_refskip() {
                    if let Some(qpos) = alignment.qpos() {
                        let a = record.seq()[qpos];
                        let q = record.qual()[qpos];

                        bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                        bmap.get_mut(&seq_name).unwrap().1.push(q + 33);
                    }
                }

                if let bam::pileup::Indel::Ins(len) = alignment.indel() {
                    if let Some(pos1) = alignment.qpos() {
                        let pos2 = pos1 + (len as usize);
                        for pos in pos1..pos2 {
                            let a = record.seq()[pos];
                            let q = record.qual()[pos];

                            bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                            bmap.get_mut(&seq_name).unwrap().1.push(q + 33);
                        }
                    }
                }
            }
        }
    }
    
    pb.finish_with_message("Pileup processing completed");

    // Filter reads that don't span the full locus (at least 80% of the region)
    let filtered_bmap: HashMap<String, (String, Vec<u8>)> = bmap.clone()
        .into_iter()
        .filter(|(seq_name, _)| {
            // Extract read name from sequence name
            if let Some(read_name) = seq_name.split('|').next() {
                if let Some((read_start, read_end)) = read_spans.get(read_name) {
                    *read_start <= start && *read_end >= end - 1
                    
                } else {
                    false
                }
            } else {
                false
            }
        })
        .collect();

    let records = filtered_bmap
        .iter()
        .map(|sq| fastq::Record::with_attrs(sq.0.as_str(), None, sq.1.0.as_bytes(), &sq.1.1))
        .collect();
    let filtered_read_span: Vec<(String, u64, u64)> = filtered_bmap.iter().map(|(seq_name, _)| {
        let read_name = seq_name.split('|').next().unwrap();
        let read_span = read_spans.get(read_name).unwrap();
        (read_name.to_string(), read_span.0, read_span.1)
    }).collect();
    // println!("read_span: {:?}", filtered_read_span);
    println!("original_bmap{:?}, filtered_bmap: {:?}", bmap.len(), filtered_bmap.len());

    Ok((records))
}
// Function to extract seqs from a FASTA file within a specified genomic region.
pub fn extract_fasta_seqs(
    basename: &String,
    fasta: &mut Reader,
    chr: &String,
    start: &u64,
    stop: &u64,
    name: &String,
) -> AnyhowResult<Vec<fastq::Record>> {
    let id = format!("{chr}:{start}-{stop}|{name}|{basename}");
    let seq = fasta.fetch_seq_string(chr, usize::try_from(*start)?, usize::try_from(*stop - 1).unwrap()).unwrap();

    if seq.len() > 0 {
        let records = vec![fastq::Record::with_attrs(id.as_str(), None, seq.as_bytes(), vec![30; seq.len()].as_slice())];

        return Ok(records);
    }

    Err(anyhow::anyhow!("No sequence found for locus: {}", id))
}

