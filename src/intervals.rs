use log::{info, warn, debug};
use anyhow::{Result as AnyhowResult, Context};
use bio::io::fastq;
use rust_htslib::faidx::Reader;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, Record, record::Aux, pileup::Pileup};
use std::path::{PathBuf};
use std::collections::{HashMap, BTreeMap, HashSet};
use std::io::Write;
use std::fs::File;
use crate::util;
use indicatif::{ProgressBar, ProgressStyle};
use crate::methyl;

pub fn get_read_start_end(pileup: &mut Pileup, primary_only: bool) -> AnyhowResult<(HashMap<String, u64>, HashMap<String, Vec<u8>>, HashMap<String, Vec<u8>>)> {
    let mut read_coordinates: HashMap<String, u64> = HashMap::new();
    let mut read_sequence_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut readnames = HashSet::new();
    let mut read_quality_dict: HashMap<String, Vec<u8>> = HashMap::new();
    for (i, alignment) in pileup.alignments().enumerate() {
        let record = alignment.record();
        
        let qname = String::from_utf8_lossy(record.qname()).into_owned() ;
        let read_seq = record.seq().clone().as_bytes();
        // skip the alignment from the same read
        if readnames.contains(&qname) {
            continue;
        }
        readnames.insert(qname.clone());

        let is_secondary = record.is_secondary();
        let is_supplementary = record.is_supplementary();
        if primary_only && (is_secondary || is_supplementary) {
            continue;
        }

        if !read_sequence_dict.contains_key(&qname) {
            read_sequence_dict.insert(qname.clone(), read_seq);
        }
        if !read_quality_dict.contains_key(&qname) {
            read_quality_dict.insert(qname.clone(), record.qual().to_vec());
        }
        
        // Handle different alignment types
        match alignment.indel() {
            bam::pileup::Indel::Ins(len) => {
                
                // For insertions, add one reference base followed by the insertion bases
                if let Some(pos1) = alignment.qpos() {
                    read_coordinates.insert(qname.clone(), pos1 as u64);
                }
            }
            bam::pileup::Indel::Del(_) => {
                // For deletions, add the first base of the deletion
                if let Some(qpos) = alignment.qpos() {
                    read_coordinates.insert(qname.clone(), qpos as u64);
                }                       
            }
            bam::pileup::Indel::None => {
                // For matches/mismatches, add the base
                if let Some(qpos) = alignment.qpos() {
                    read_coordinates.insert(qname.clone(), qpos as u64);
                }
            
            }
        }
    }
    Ok((read_coordinates, read_sequence_dict, read_quality_dict))
}


// extract haplotypes from bam file
pub fn extract_haplotypes_from_bam(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    sampleid: &String,
    primary_only: bool,
) -> AnyhowResult<(Vec<fastq::Record>)> {
    let rg_sm_map = util::get_rg_to_sm_mapping(bam);
    let mut read_start: HashMap<String, u64> = HashMap::new();
    let mut read_end: HashMap<String, u64> = HashMap::new();
    let mut read_sequence_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut read_quality_dict: HashMap<String, Vec<u8>> = HashMap::new();
    // Create progress bar
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap()
    );
    pb.set_message("Processing pileup...");

    let _ = bam.fetch((chr.as_bytes(), start, end));

    for p in bam.pileup() {
        let mut pileup = p?;
        if pileup.pos() as u64 == start {
            let (read_s, read_sequence_s_d, read_quality_s_d) = get_read_start_end(&mut pileup, primary_only)?;
            read_start.extend(read_s);
            read_sequence_dict.extend(read_sequence_s_d);
            read_quality_dict.extend(read_quality_s_d);
        }
        if pileup.pos() as u64 == end {
            let (read_e, read_sequence_e_d, read_quality_e_d) = get_read_start_end(&mut pileup, primary_only)?;
            read_end.extend(read_e);
            read_sequence_dict.extend(read_sequence_e_d);
            read_quality_dict.extend(read_quality_e_d);
        }
    }
    
    let mut filtered_sequence_dict = HashMap::new();
    for (read_name, seq) in read_sequence_dict.iter() {
        if read_start.contains_key(read_name) && read_end.contains_key(read_name) {
            let read_start = read_start.get(read_name).unwrap();
            let read_end = read_end.get(read_name).unwrap();
            let read_seq = seq.clone();
            let read_seq_sub = read_seq[*read_start as usize..*read_end as usize + 1].to_vec();
            let read_quality = read_quality_dict.get(read_name).unwrap();
            let read_quality_sub = read_quality[*read_start as usize..*read_end as usize + 1].to_vec();
            let record_id = format!("{read_name}|{chr}:{start}-{end}|{sampleid}");
            filtered_sequence_dict.insert(record_id.clone(), (read_seq_sub, read_quality_sub));
        }
    }

    let records: Vec<fastq::Record> = filtered_sequence_dict
        .iter()
        .map(|(read_name, (read_seq, read_qual))| fastq::Record::with_attrs(read_name, None, read_seq, read_qual))
        .collect();

    debug!("Extracted {} reads from region", records.len());
   
    Ok((records))
}


// extract haplotypes from bam file
pub fn extract_haplotypes_coordinates_from_bam(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    primary_only: bool,
) -> AnyhowResult<(Vec<fastq::Record>, HashMap<String, (u64, u64)>, HashMap<String, String>, HashMap<String, Record>)> {
    let rg_sm_map = util::get_rg_to_sm_mapping(bam);
    let mut bmap: HashMap<String, (String, Vec<u8>)> = HashMap::new();
    let mut read_spans: HashMap<String, (u64, u64)> = HashMap::new(); // Track read alignment spans
    let mut read_coordinates: HashMap<String, (u64, u64)> = HashMap::new(); // Track read coordinates
    let _ = bam.fetch((chr.as_bytes(), start, end));
    let mut read_sequence_dict: HashMap<String, String> = HashMap::new();
    let mut bam_records: HashMap<String, Record> = HashMap::new();
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
                bam_records.insert(qname.clone(), record.clone());

                let sm = match util::get_sm_name_from_rg(&record, &rg_sm_map) {
                    Ok(a) => a,
                    Err(_) => String::from("unknown"),
                };

                let is_secondary = record.is_secondary();
                let is_supplementary = record.is_supplementary();
                if primary_only && (is_secondary || is_supplementary) {
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
                        
                        // For insertions, add one reference base followed by the insertion bases
                        if let Some(pos1) = alignment.qpos() {
                            // Then add the insertion bases
                            let pos2 = pos1 + (len as usize) + 1;
                            for pos in pos1..pos2 {
                                let a = record.seq()[pos];
                                let q = record.qual()[pos];
                                let valid_q = q.min(40);

                                // Track coordinates in the reconstructed sequence
                                let read_start = read_coordinates.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
                                read_start.0 = read_start.0.min(pos as u64);
                                read_start.1 = read_start.1.max(pos as u64);

                                bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                                bmap.get_mut(&seq_name).unwrap().1.push(valid_q + 33);
                            }
                        }
                    }
                    bam::pileup::Indel::Del(_) => {
                        // For deletions, add the first base of the deletion
                        if let Some(qpos) = alignment.qpos() {
                            let a = record.seq()[qpos];
                            let q = record.qual()[qpos];
                            let valid_q = q.min(40);
                            
                            
                            // Track coordinates in the reconstructed sequence
                            let read_start = read_coordinates.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
                            read_start.0 = read_start.0.min(qpos as u64);
                            read_start.1 = read_start.1.max(qpos as u64);

                            bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                            bmap.get_mut(&seq_name).unwrap().1.push(valid_q+33);
                        }                       
                    }
                    bam::pileup::Indel::None => {
                        // For matches/mismatches, add the base

                        if let Some(qpos) = alignment.qpos() {
                            let a = record.seq()[qpos];
                            let q = record.qual()[qpos];
                            let valid_q = q.min(40);
                            
                            // Track coordinates in the reconstructed sequence
                            let read_start = read_coordinates.entry(qname.clone()).or_insert((u64::MAX, u64::MIN));
                            read_start.0 = read_start.0.min(qpos as u64);
                            read_start.1 = read_start.1.max(qpos as u64);

                            bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                            bmap.get_mut(&seq_name).unwrap().1.push(valid_q + 33);
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
    // println!("filtered_bmap: {}", filtered_bmap.len());

    let records: Vec<fastq::Record> = filtered_bmap
        .iter()
        .map(|sq| fastq::Record::with_attrs(sq.0.as_str(), None, sq.1.0.as_bytes(), sq.1.1.as_slice()))
        .collect();

    debug!("Extracted {} reads from region", records.len());

   
    Ok((records, read_coordinates, read_sequence_dict, bam_records))
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
    if !seq.is_empty() {
        let records = vec![fastq::Record::with_attrs(id.as_str(), None, seq.as_bytes(), vec![30; seq.len()].as_slice())];

        return Ok(records);
    }

    Err(anyhow::anyhow!("No sequence found for locus: {}", id))
}

pub fn process_bam_file_by_coordinates(bam: &mut IndexedReader, chromosome: &str, start: usize, end: usize, primary_only: bool, sampleid: &String) -> Vec<fastq::Record> {
    
    let reads = extract_haplotypes_from_bam(
        bam,
        &chromosome,
        start as u64,
        end as u64,
        &sampleid,
        primary_only, 
    ).unwrap();
    debug!("Extracted {} reads from region", reads.len());
    reads
}


pub fn process_bam_file(bam: &mut IndexedReader, chromosome: &str, start: usize, end: usize, primary_only: bool) -> (Vec<fastq::Record>, HashMap<String,HashMap<usize, f32>>) {
    
    let (reads, read_coordinates, read_sequence_dictionary, bam_records) = extract_haplotypes_coordinates_from_bam(
        bam,
        &chromosome,
        start as u64,
        end as u64,
        primary_only, 
    ).unwrap();
    let methyl_all_reads = methyl::start(bam_records,  &read_coordinates);
    debug!("Extracted {} reads from region", reads.len());
    let mut filtered_reads = Vec::new();
    for r in reads.iter() {
        let r_name = r.id().to_string().split('|').next().unwrap().to_string();
        let r_seq = read_sequence_dictionary.get(&r_name).unwrap();
        if read_coordinates.contains_key(&r_name){  
            let r_start = read_coordinates.get(&r_name).unwrap().0;
            let r_end = read_coordinates.get(&r_name).unwrap().1 + 1;
            let reconstructed_seq = String::from_utf8_lossy(r.seq()).to_string();
            // println!("r_seq: {}, reconstructed_seq: {}, readname: {}", r_seq.contains(reconstructed_seq.as_str()), reconstructed_seq, r_name );
            let r_sub_seq = r_seq[r_start as usize..r_end as usize].to_string();
            // sanity check if the reconstructed sequence is the same as the read sequence
            if r_sub_seq == reconstructed_seq{
                filtered_reads.push(r.clone())
            }else{
                println!("r_seq: {}, reconstructed_seq: {}, readname: {}", r_sub_seq, reconstructed_seq, r_name );
            }
        }else{
            // large deletion in the whole region
            assert!(r.seq().len() < 1);
            filtered_reads.push(r.clone())
        }
        
    }
    
    (filtered_reads, methyl_all_reads)
}

pub fn process_fasta_file(reference: &Vec<fastq::Record>, chromosome: &str, start:usize, end:usize, sampleid: &String) -> Vec<fastq::Record> {
    let mut reference_seqs = Vec::new();
    for record in reference.iter() {
        let seq_id = record.id().to_string();
        let sequence = String::from_utf8_lossy(record.seq()).to_string();
        // Check if this sequence matches the target chromosome
        if seq_id == chromosome {
            debug!("Extracting region {}:{}-{} from chromosome {}", 
                chromosome, start, end, seq_id);
            // Extract the specified region
            if start < sequence.len() && end <= sequence.len() {
                let region_seq = sequence[start..end].to_string();
                let record_id = format!("{}:{}-{}|reference|{}", 
                                    chromosome, start, end, sampleid);
                
                let fastq_record = fastq::Record::with_attrs(
                    &record_id,
                    None,
                    region_seq.as_bytes(),
                    vec![30; region_seq.len()].as_slice() // Default quality score
                );
                
                reference_seqs.push(fastq_record);
                debug!("Extracted reference sequence: {} bp", region_seq.len());
            } else {
                warn!("Region coordinates out of bounds for sequence length {}", sequence.len());
            }
            
        }
    }
    
    if reference_seqs.is_empty() {
        warn!("No matching chromosome '{}' found in FASTA file", chromosome);
    }
    
    debug!("Extracted {} reference sequences", reference_seqs.len());
    reference_seqs
}


pub fn write_fasta_output(final_hap: HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)>, output_path: &PathBuf) -> AnyhowResult<()> {

    
    let mut file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;
    
    let mut index:i32 = 1;
    let chars_per_line = 60;

    for (sequence, (cigar, read_names, allele_frequency)) in final_hap {
        // Use the first read name as the sequence ID, or create a hash-based ID
        let sequence_id = if !read_names.is_empty() {
            let json_string =serde_json::to_string(&read_names).expect("Failed to serialize to JSON");
            json_string
        } else {
            "unknown".to_string()   
        };
        writeln!(file, ">Haplotype {}\tLength{}\tCigar {}\tSupportReads {}\tAlleleFrequency {:.2}\t{}", index, sequence.len(), cigar, read_names.len(), allele_frequency, sequence_id)?;
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

        index += 1;
    }


    
    Ok(())
}

pub fn collapse_haplotypes(reads: &Vec<fastq::Record>, read_methyl_dict: &HashMap<String,HashMap<usize, f32>>, reference: &fastq::Record, min_reads: usize, frequency_min: f64) -> AnyhowResult<HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)>> {
    // Read the reads records (name and sequence) into a vector.
    let mut read_dictionary: HashMap<String, HashMap<String, HashMap<usize, f32>>> = HashMap::new();
    for read in reads{
        let read_name_origin = read.id().to_string();
        let contig_name = read_name_origin.split('|').next().unwrap().to_string();
        // println!("contig_name: {:?}", contig_name);
        let contig = String::from_utf8_lossy(read.seq()).to_string();
        if !read_methyl_dict.contains_key(&contig_name){
            read_dictionary.entry(contig).or_default().insert(contig_name, HashMap::new());
        }else{
            let methyl_dict = read_methyl_dict.get(&contig_name).unwrap();
            read_dictionary.entry(contig).or_default().insert(contig_name, methyl_dict.clone());
        }
    }

    let reference_seq = String::from_utf8_lossy(reference.seq()).to_string();
    let mut final_hap: HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)> = HashMap::new();

    for (hap, vec) in read_dictionary.iter() {
        let allele_frequency = vec.len() as f64 / reads.len() as f64;
        if (vec.len() >= min_reads as usize) && (allele_frequency >= frequency_min as f64){
            let cigar = util::gap_open_aligner(&reference_seq, &hap);
            final_hap.insert(hap.clone(), (cigar, vec.clone(), allele_frequency));
        }
    }
    debug!("Extracted {} haplotypes from region, total reads: {:?}", final_hap.len(), read_dictionary.values().map(|x| x.len()).collect::<Vec<_>>());
    Ok(final_hap)
}




pub fn start(bam: &mut IndexedReader, reference_fa: &Vec<fastq::Record>, chromosome: &str, start: usize, end: usize, sampleid: &String, min_reads: usize, frequency_min: f64, primary_only: bool, write_output:bool) -> AnyhowResult<HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)>>  {
    let (reads_list, read_methyl_dict) = process_bam_file( bam, chromosome, start, end, primary_only);  
    // let reads_list = process_bam_file_by_coordinates(bam, chromosome, start, end, primary_only, sampleid);
    let reference = process_fasta_file( reference_fa, chromosome, start, end, sampleid);
    let reference = reference.first().unwrap().clone();
    let final_hap = collapse_haplotypes(&reads_list, &read_methyl_dict, &reference, min_reads, frequency_min)?;
    info!("Haplotype reconstruction completed");
    let final_hap_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplograph.fasta", ".", sampleid, chromosome, start, end));
    if write_output{    
        write_fasta_output(final_hap.clone(), &final_hap_output)?;
    }
   
    Ok((final_hap))
}