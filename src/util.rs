use std::collections::{HashMap, HashSet};
use bio::io::fastq;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use rust_htslib::faidx::Reader;
use anyhow::{Result as AnyhowResult};
use indicatif::{ProgressBar, ProgressStyle};
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;

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



pub fn split_locus(locus: String) -> (String, usize, usize) {
    let parts: Vec<&str> = locus.split(':').collect();
    let chromosome = parts[0].to_string();
    let start: usize = parts[1].split("-").collect::<Vec<&str>>().first().unwrap().parse().unwrap();
    let end: usize = parts[1].split("-").collect::<Vec<&str>>().last().unwrap().parse().unwrap();
    (chromosome, start, end)
}

pub fn extract_haplotypes_coordinates_from_bam(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    min_mapping_quality: Option<u8>,
) -> AnyhowResult<(Vec<fastq::Record>, HashMap<String, (u64, u64)>, HashMap<String, String>)> {
    let rg_sm_map = get_rg_to_sm_mapping(bam);
    let mut bmap = HashMap::new();
    let mut read_spans: HashMap<String, (u64, u64)> = HashMap::new(); // Track read alignment spans
    let mut read_coordinates: HashMap<String, (u64, u64)> = HashMap::new(); // Track read coordinates
    let _ = bam.fetch((chr.as_bytes(), start, end));
    let mut read_sequence_dict: HashMap<String, String> = HashMap::new();
    // Create progress bar
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap()
    );
    pb.set_message("Processing pileup...");
    // let mut insertion_positions = HashMap::new();
    for p in bam.pileup() {
        let pileup = p?;
        pb.set_message(format!("Processing position {}", pileup.pos()));
        let mut readnames = HashSet::new();
        if start <= (pileup.pos() as u64) && (pileup.pos() as u64) < end {
            for (i, alignment) in pileup.alignments().enumerate() {
                let record = alignment.record();
                
                let qname = String::from_utf8_lossy(record.qname()).into_owned() ;
                let read_seq = String::from_utf8_lossy(&record.seq().as_bytes()).into_owned();
                // skip the alignment from the same read
                if readnames.contains(&qname) {
                    continue;
                }
                readnames.insert(qname.clone());

                let sm = match get_sm_name_from_rg(&record, &rg_sm_map) {
                    Ok(a) => a,
                    Err(_) => String::from("unknown"),
                };

                let is_secondary = record.is_secondary();
                let is_supplementary = record.is_supplementary();
                if is_secondary || is_supplementary {
                    continue;
                }
                if !read_sequence_dict.contains_key(&qname) {
                    read_sequence_dict.insert(qname.clone(), read_seq);
                }
                

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

                // Handle different alignment types
                match alignment.indel() {
                    bam::pileup::Indel::Ins(len) => {
                        // insertion_positions.get_mut(&seq_name).unwrap().push((pileup.pos() as u64, len as u64));
                        
                        // For insertions, add one reference base followed by the insertion bases
                        if let Some(pos1) = alignment.qpos() {
                            // Then add the insertion bases
                            let pos2 = pos1 + (len as usize) + 1;
                            for pos in pos1..pos2 {
                                let a = record.seq()[pos];
                                let q = record.qual()[pos];

                                // Track coordinates in the reconstructed sequence
                                let read_start = read_coordinates.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
                                read_start.0 = read_start.0.min(pos as u64);
                                read_start.1 = read_start.1.max(pos as u64);

                                bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                                bmap.get_mut(&seq_name).unwrap().1.push(q + 33);
                            }
                        }
                    }
                    bam::pileup::Indel::Del(_) => {
                        // For deletions, add the first base of the deletion
                        if let Some(qpos) = alignment.qpos() {
                            let a = record.seq()[qpos];
                            let q = record.qual()[qpos];
                            
                            
                            // Track coordinates in the reconstructed sequence
                            let read_start = read_coordinates.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
                            read_start.0 = read_start.0.min(qpos as u64);
                            read_start.1 = read_start.1.max(qpos as u64);

                            bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                            bmap.get_mut(&seq_name).unwrap().1.push(q + 33);
                        }                       
                    }
                    bam::pileup::Indel::None => {
                        // For matches/mismatches, add the base

                        if let Some(qpos) = alignment.qpos() {
                            let a = record.seq()[qpos];
                            let q = record.qual()[qpos];
                            
                            // Track coordinates in the reconstructed sequence
                            let read_start = read_coordinates.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
                            read_start.0 = read_start.0.min(qpos as u64);
                            read_start.1 = read_start.1.max(qpos as u64);

                            bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                            bmap.get_mut(&seq_name).unwrap().1.push(q + 33);
                        }
                    
                    }
                }
            }
        }
    }
    // println!("insertion_positions: {:?}", insertion_positions);
    
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

    let records: Vec<fastq::Record> = filtered_bmap
        .iter()
        .map(|sq| fastq::Record::with_attrs(sq.0.as_str(), None, sq.1.0.as_bytes(), &sq.1.1))
        .collect();


   
    Ok((records, read_coordinates, read_sequence_dict))
}

// pub fn extract_haplotypes_from_bam(
//     bam: &mut IndexedReader,
//     chr: &str,
//     start: u64,
//     end: u64,
//     min_mapping_quality: Option<u8>,
// ) -> AnyhowResult<(Vec<fastq::Record>, HashMap<String, (u64, u64)>, HashMap<String, String>)> {
//     let rg_sm_map = get_rg_to_sm_mapping(bam);
//     let mut bmap = HashMap::new();
//     let mut read_spans: HashMap<String, (u64, u64)> = HashMap::new(); // Track read alignment spans
//     let mut read_coordinates: HashMap<String, (u64, u64)> = HashMap::new(); // Track read coordinates
//     let _ = bam.fetch((chr.as_bytes(), start, end));
//     let mut read_sequence_dict: HashMap<String, String> = HashMap::new();
//     // Create progress bar
//     let pb = ProgressBar::new_spinner();
//     pb.set_style(
//         ProgressStyle::default_spinner()
//             .template("{spinner:.green} [{elapsed_precise}] {msg}")
//             .unwrap()
//     );
//     pb.set_message("Processing pileup...");
    
//     for p in bam.pileup() {
//         let pileup = p?;
//         pb.set_message(format!("Processing position {}", pileup.pos()));

//         if start <= (pileup.pos() as u64) && (pileup.pos() as u64) < end {
//             for (i, alignment) in pileup.alignments().enumerate() {
//                 let record = alignment.record();
                
//                 let qname = String::from_utf8_lossy(record.qname()).into_owned() ;
//                 let read_seq = String::from_utf8_lossy(&record.seq().as_bytes()).into_owned();

//                 let sm = match get_sm_name_from_rg(&record, &rg_sm_map) {
//                     Ok(a) => a,
//                     Err(_) => String::from("unknown"),
//                 };

//                 let is_secondary = record.is_secondary();
//                 let is_supplementary = record.is_supplementary();
//                 if is_secondary || is_supplementary {
//                     continue;
//                 }
//                 if !read_sequence_dict.contains_key(&qname) {
//                     read_sequence_dict.insert(qname.clone(), read_seq);
//                 }
                

//                 let locus_key = format!("{}:{}-{}", chr, start, end);
//                 let seq_name = format!("{qname}|{locus_key}|{sm}");
//                 // let seq_name = qname.clone();

//                 // Track alignment span for this read
//                 let reference_start = record.pos() as u64;
//                 assert!(record.cigar().pos() == reference_start as i64);
//                 let reference_end = record.cigar().end_pos() as u64; 
//                 let read_span = read_spans.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
//                 read_span.0 = read_span.0.min(reference_start);
//                 read_span.1 = read_span.1.max(reference_end);

//                 if !bmap.contains_key(&seq_name) {
//                     bmap.insert(seq_name.clone(), (String::new(), Vec::new()));
//                 }

//                 // Handle different alignment types
//                 match alignment.indel() {
//                     bam::pileup::Indel::Ins(len) => {
//                         // For insertions, add one reference base followed by the insertion bases
//                         if let Some(pos1) = alignment.qpos() {
//                             // Then add the insertion bases
//                             let pos2 = pos1 + (len as usize) + 1;
//                             for pos in pos1..pos2 {
//                                 let a = record.seq()[pos];
//                                 let q = record.qual()[pos];

//                                 // Track coordinates in the reconstructed sequence
//                                 let read_start = read_coordinates.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
//                                 read_start.0 = read_start.0.min(pos as u64);
//                                 read_start.1 = read_start.1.max(pos as u64);

//                                 bmap.get_mut(&seq_name).unwrap().0.push(a as char);
//                                 bmap.get_mut(&seq_name).unwrap().1.push(q + 33);
//                             }
//                         }
//                     }
//                     bam::pileup::Indel::Del(_) => {
//                         // For deletions, add the first base of the deletion
//                         if let Some(qpos) = alignment.qpos() {
//                             let a = record.seq()[qpos];
//                             let q = record.qual()[qpos];
                            
                            
//                             // Track coordinates in the reconstructed sequence
//                             let read_start = read_coordinates.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
//                             read_start.0 = read_start.0.min(qpos as u64);
//                             read_start.1 = read_start.1.max(qpos as u64);

//                             bmap.get_mut(&seq_name).unwrap().0.push(a as char);
//                             bmap.get_mut(&seq_name).unwrap().1.push(q + 33);
//                         }                       
//                     }
//                     bam::pileup::Indel::None => {
//                         // For matches/mismatches, add the base

//                         if let Some(qpos) = alignment.qpos() {
//                             let a = record.seq()[qpos];
//                             let q = record.qual()[qpos];
                            
//                             // Track coordinates in the reconstructed sequence
//                             let read_start = read_coordinates.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
//                             read_start.0 = read_start.0.min(qpos as u64);
//                             read_start.1 = read_start.1.max(qpos as u64);

//                             bmap.get_mut(&seq_name).unwrap().0.push(a as char);
//                             bmap.get_mut(&seq_name).unwrap().1.push(q + 33);
//                         }
                    
//                     }
//                 }
//             }
//         }
//     }
    
//     pb.finish_with_message("Pileup processing completed");

//     // Filter reads that don't span the full locus (at least 80% of the region)
//     let filtered_bmap: HashMap<String, (String, Vec<u8>)> = bmap.clone()
//         .into_iter()
//         .filter(|(seq_name, _)| {
//             // Extract read name from sequence name
//             if let Some(read_name) = seq_name.split('|').next() {
//                 if let Some((read_start, read_end)) = read_spans.get(read_name) {
//                     *read_start <= start && *read_end >= end - 1
                    
//                 } else {
//                     false
//                 }
//             } else {
//                 false
//             }
//         })
//         .collect();

//     let records: Vec<fastq::Record> = filtered_bmap
//         .iter()
//         .map(|sq| fastq::Record::with_attrs(sq.0.as_str(), None, sq.1.0.as_bytes(), &sq.1.1))
//         .collect();


   
//     Ok((records, read_coordinates, read_sequence_dict))
// }
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

    if !seq.is_empty() {
        let records = vec![fastq::Record::with_attrs(id.as_str(), None, seq.as_bytes(), vec![30; seq.len()].as_slice())];

        return Ok(records);
    }

    Err(anyhow::anyhow!("No sequence found for locus: {}", id))
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


// /// Extract read names and sequences from a specific region in a BAM file
// pub fn extract_read_sequences_from_region(
//     bam: &mut IndexedReader,
//     chr: &str,
//     start: u64,
//     end: u64,
// ) -> AnyhowResult<HashMap<String, String>> {
//     let mut read_sequences: HashMap<String, String> = HashMap::new();
    
//     // Fetch the specific region
//     let _ = bam.fetch((chr.as_bytes(), start, end));
    
//     for result in bam.records() {
//         let record = result?;
//         let read_name = String::from_utf8_lossy(record.qname()).into_owned();
//         let read_seq = String::from_utf8_lossy(&record.seq().as_bytes()).into_owned();
//         if read_sequences.contains_key(&read_name) {
//             // println!("read_name: {} already exists", read_name);
//             // continue;
//             let read_seq_indict = read_sequences.get(&read_name).unwrap().clone();
//             if read_seq !=  read_seq_indict {
//                 println!("read_name: {} already exists, but sequences are different, {} != {}", read_name, read_seq.len(), read_seq_indict.len());
//             }
//         }

        
//         read_sequences.insert(read_name, read_seq);
    
//     }
    
//     Ok(read_sequences)
// }

/// Find overlapping reads between two read vectors
pub fn find_overlapping_reads(read_vector1: &[String], read_vector2: &[String]) -> Vec<String> {
    let set1: std::collections::HashSet<String> = read_vector1.iter().map(|s| s.split('|').next().unwrap().to_string()).collect();
    let set2: std::collections::HashSet<String> = read_vector2.iter().map(|s| s.split('|').next().unwrap().to_string()).collect();
    
    set1.intersection(&set2)
        .map(|read_name| read_name.clone())
        .collect()
}

