use crate::intervals;
use anyhow::{Context, Result as AnyhowResult};
use rust_htslib::bam::IndexedReader;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::collections::HashMap;
use std::collections::HashSet;
use rust_htslib::bam::{self, Read as BamRead, pileup::Pileup};
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use log::info;
use crate::util;


pub fn get_read_start_end(pileup: &mut Pileup, primary_only: bool) -> AnyhowResult<(HashMap<String, u64>, HashMap<String, Vec<u8>>, HashMap<String, Vec<u8>>, HashMap<String, String>)> {
    let mut read_coordinates: HashMap<String, u64> = HashMap::new();
    let mut read_sequence_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut readnames = HashSet::new();
    let mut read_quality_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut read_strand_dict: HashMap<String, String> = HashMap::new();
    for (i, alignment) in pileup.alignments().enumerate() {
        let record = alignment.record();

        let qname = String::from_utf8_lossy(record.qname()).into_owned() ;
        let read_seq = record.seq().clone().as_bytes();
        // skip the alignment from the same read
        if readnames.contains(&qname) {
            continue;
        }
        readnames.insert(qname.clone());

        let strand = record.strand().to_string();
        if !read_strand_dict.contains_key(&qname) {
            read_strand_dict.insert(qname.clone(), strand);
        }

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
    Ok((read_coordinates, read_sequence_dict, read_quality_dict, read_strand_dict))
}


// extract haplotypes from bam file
pub fn extract_haplotypes_from_bam(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    sampleid: &String,
    primary_only: bool,
) -> AnyhowResult<HashMap<String, String>> {
    let mut read_start: HashMap<String, u64> = HashMap::new();
    let mut read_end: HashMap<String, u64> = HashMap::new();
    let mut read_sequence_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut read_quality_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut read_strand_dict: HashMap<String, String> = HashMap::new();
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
            let (read_s, read_sequence_s_d, read_quality_s_d, read_strand_s_d) = get_read_start_end(&mut pileup, primary_only)?;
            read_start.extend(read_s);
            read_sequence_dict.extend(read_sequence_s_d);
            read_quality_dict.extend(read_quality_s_d);
            read_strand_dict.extend(read_strand_s_d);
        }
        if pileup.pos() as u64 == end {
            let (read_e, read_sequence_e_d, read_quality_e_d, read_strand_e_d) = get_read_start_end(&mut pileup, primary_only)?;
            read_end.extend(read_e);
            read_sequence_dict.extend(read_sequence_e_d);
            read_quality_dict.extend(read_quality_e_d);
            read_strand_dict.extend(read_strand_e_d);
        }
    }
    
    let mut filtered_sequence_dict = HashMap::new();
    for (read_name, seq) in read_sequence_dict.iter() {
        if read_start.contains_key(read_name) && read_end.contains_key(read_name) {
            let read_start = *read_start.get(read_name).unwrap() as usize;
            let read_end = *read_end.get(read_name).unwrap() as usize;
            let read_seq = seq.clone();
            let pos_start = read_start.min(read_end);
            let pos_end = read_start.max(read_end);
            let read_strand = read_strand_dict.get(read_name).unwrap();

            let read_seq_sub = String::from_utf8_lossy(&read_seq[pos_start..pos_end]).to_string();
            let read_seq_final = if read_strand == "+" {
                read_seq_sub
            } else {
                util::reverse_complement(&read_seq_sub)
            };

            let record_id = format!("{read_name}|{chr}:{start}-{end}|{sampleid}");
            filtered_sequence_dict.insert(record_id.clone(), read_seq_final);
        }
    }

    info!("Extracted {} reads from region", filtered_sequence_dict.len());
   
    Ok(filtered_sequence_dict)
}


pub fn start(
    bam: &mut IndexedReader,
    chromosome: &str,
    start: usize,
    end: usize,
    primary_only: bool,
    output_path: String,
    sampleid: String,
    pileup:bool,
) -> AnyhowResult<()> {
    if pileup{
        let (reads, read_coordinates, read_sequence_dictionary, bam_records) =
        intervals::extract_haplotypes_coordinates_from_bam(
            bam,
            chromosome,
            start as u64,
            end as u64,
            primary_only,
        )
        .unwrap();
        let outputfile = PathBuf::from(format!("{}.fasta", output_path));
        let mut file = File::create(outputfile)
            .with_context(|| format!("Failed to create output file: {}.fasta", output_path))?;

        for r in reads.iter() {
            let r_name = String::from_utf8_lossy(r.id().as_bytes()).to_string();
            let r_seq = String::from_utf8_lossy(r.seq()).to_string();
            writeln!(file, ">{}", r_name)?;
            // write the sequence in fasta format
            let seq_len = r_seq.len();
            let chars_per_line = 60;
            let full_lines = seq_len / chars_per_line;
            for i in 0..full_lines {
                let start = i * chars_per_line;
                let end = start + chars_per_line;
                writeln!(file, "{}", &r_seq[start..end])?;
            }
            if seq_len % chars_per_line != 0 {
                writeln!(file, "{}", &r_seq[full_lines * chars_per_line..])?;
            }
        }
    }else {
        let reads = extract_haplotypes_from_bam(
            bam,
            &chromosome,
            start as u64,
            end as u64,
            &sampleid,
            primary_only, 
        ).unwrap();

        let output_p = PathBuf::from(format!("{}.fasta", output_path));
        if ! reads.is_empty(){
            let _ = util::write_fasta(&reads, &output_p);
        }  
    }


    Ok(())
}
