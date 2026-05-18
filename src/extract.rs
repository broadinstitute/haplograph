use crate::intervals;
use anyhow::{Context, Result as AnyhowResult};
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, pileup::Pileup, Record as BamRecord};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::collections::HashMap;
use std::collections::HashSet;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use log::info;
use crate::util;


pub fn get_read_start_end(pileup: &mut Pileup, primary_only: bool) -> AnyhowResult<(HashMap<String, u64>, HashMap<String, Vec<u8>>, HashMap<String, Vec<u8>>, HashMap<String, String>, HashMap<String, BamRecord>)> {
    let mut read_coordinates: HashMap<String, u64> = HashMap::new();
    let mut read_sequence_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut readnames = HashSet::new();
    let mut read_quality_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut read_strand_dict: HashMap<String, String> = HashMap::new();
    let mut bam_records_dict: HashMap<String, BamRecord> = HashMap::new();
    for (i, alignment) in pileup.alignments().enumerate() {
        let record = alignment.record();

        let qname = String::from_utf8_lossy(record.qname()).into_owned() ;
        bam_records_dict.insert(qname.clone(), record.clone());
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
    Ok((read_coordinates, read_sequence_dict, read_quality_dict, read_strand_dict, bam_records_dict))
}

pub fn normalized_read_slice_bounds(
    read_start: usize,
    read_end: usize,
    read_len: usize,
    read_strand: &str,
) -> Option<(usize, usize)> {
    let (pos_start, pos_end) = if read_strand == "+" {
        (read_start.min(read_end), read_start.max(read_end))
    } else {
        (
            read_len.saturating_sub(read_start.max(read_end)),
            read_len.saturating_sub(read_start.min(read_end)),
        )
    };
    if pos_end <= pos_start || pos_end > read_len {
        return None;
    }
    Some((pos_start, pos_end))
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

    let _ = bam.fetch((chr.as_bytes(), start-5, end+5));

    for p in bam.pileup() {
        let mut pileup = p?;
        if pileup.pos() as u64 == start {
            let (read_s, read_sequence_s_d, read_quality_s_d, read_strand_s_d, _) = get_read_start_end(&mut pileup, primary_only)?;
            read_start.extend(read_s);
            read_sequence_dict.extend(read_sequence_s_d);
            read_quality_dict.extend(read_quality_s_d);
            read_strand_dict.extend(read_strand_s_d);
        }
        if pileup.pos() as u64 == end {
            let (read_e, read_sequence_e_d, read_quality_e_d, read_strand_e_d, _) = get_read_start_end(&mut pileup, primary_only)?;
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
            let read_strand = read_strand_dict.get(read_name).unwrap();
            let read_len = read_seq.len();
            let Some((pos_start, pos_end)) =
                normalized_read_slice_bounds(read_start, read_end, read_len, read_strand)
            else {
                continue;
            };

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

// extract haplotypes from bam file
pub fn extract_haplotypes_coordinates_from_bam(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    sampleid: &String,
    primary_only: bool,
) -> AnyhowResult<(HashMap<String, (u64, u64)>, HashMap<String, Vec<u8>>, HashMap<String, Vec<u8>>, HashMap<String, String>, HashMap<String, BamRecord>)> {
    let mut read_start: HashMap<String, u64> = HashMap::new();
    let mut read_end: HashMap<String, u64> = HashMap::new();
    let mut read_sequence_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut read_quality_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut read_strand_dict: HashMap<String, String> = HashMap::new();
    let mut bam_records_dict: HashMap<String, BamRecord> = HashMap::new();
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
            let (read_s, read_sequence_s_d, read_quality_s_d, read_strand_s_d, bam_records_s_d) = get_read_start_end(&mut pileup, primary_only)?;
            read_start.extend(read_s);
            read_sequence_dict.extend(read_sequence_s_d);
            read_quality_dict.extend(read_quality_s_d);
            read_strand_dict.extend(read_strand_s_d);
            bam_records_dict.extend(bam_records_s_d);
        }
        if pileup.pos() as u64 == end {
            let (read_e, read_sequence_e_d, read_quality_e_d, read_strand_e_d, bam_records_e_d) = get_read_start_end(&mut pileup, primary_only)?;
            read_end.extend(read_e);
            read_sequence_dict.extend(read_sequence_e_d);
            read_quality_dict.extend(read_quality_e_d);
            read_strand_dict.extend(read_strand_e_d);
            bam_records_dict.extend(bam_records_e_d);
        }
    }
    
    let mut read_coordinates_formatted = HashMap::new();
    let mut read_sequence_dict_formatted = HashMap::new();
    let mut read_quality_dict_formatted = HashMap::new();
    let mut read_strand_dict_formatted = HashMap::new();
    let mut bam_records_dict_formatted = HashMap::new();
    for (read_name, seq) in read_sequence_dict.iter() {
        if read_start.contains_key(read_name) && read_end.contains_key(read_name) {
            let read_start = *read_start.get(read_name).unwrap() as usize;
            let read_end = *read_end.get(read_name).unwrap() as usize;
            let read_seq = seq.clone();
            let read_strand = read_strand_dict.get(read_name).unwrap();
            let read_len = read_seq.len();
            let Some((pos_start, pos_end)) =
                normalized_read_slice_bounds(read_start, read_end, read_len, read_strand)
            else {
                continue;
            };

            let read_seq_sub = String::from_utf8_lossy(&read_seq[pos_start..pos_end]).to_string();
            let read_seq_final = if read_strand == "+" {
                read_seq_sub
            } else {
                util::reverse_complement(&read_seq_sub)
            };

            let record_id = intervals::generate_read_name(read_name, chr, start, end, sampleid);
            read_coordinates_formatted.insert(record_id.clone(), (read_start as u64, read_end as u64));
            read_sequence_dict_formatted.insert(record_id.clone(), read_seq_final.into_bytes());
            if let Some(read_qual) = read_quality_dict.get(read_name) {
                read_quality_dict_formatted.insert(record_id.clone(), read_qual.clone());
            }
            read_strand_dict_formatted.insert(record_id.clone(), read_strand.clone());
            if let Some(record) = bam_records_dict.get(read_name) {
                bam_records_dict_formatted.insert(record_id.clone(), record.clone());
            }
        }
    }

    info!("Extracted {} reads from region", read_coordinates_formatted.len());
   
    Ok((
        read_coordinates_formatted,
        read_quality_dict_formatted,
        read_sequence_dict_formatted,
        read_strand_dict_formatted,
        bam_records_dict_formatted,
    ))
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

#[cfg(test)]
mod tests {
    use super::normalized_read_slice_bounds;

    #[test]
    fn forward_strand_bounds_are_ordered() {
        let (start, end) = normalized_read_slice_bounds(10, 50, 100, "+").unwrap();
        assert!(end > start);
        assert_eq!((start, end), (10, 50));
    }

    #[test]
    fn reverse_strand_bounds_are_ordered_after_conversion() {
        let (start, end) = normalized_read_slice_bounds(80, 20, 100, "-").unwrap();
        assert!(end > start);
        assert_eq!((start, end), (20, 80));
    }

    #[test]
    fn reverse_strand_equal_coordinates_returns_none() {
        let bounds = normalized_read_slice_bounds(42, 42, 100, "-");
        assert!(bounds.is_none());
    }

    #[test]
    fn all_valid_cases_have_end_larger_than_start() {
        let strands = ["+", "-"];
        let read_len = 100usize;
        for strand in strands {
            for read_start in 0..read_len {
                for read_end in 0..read_len {
                    if let Some((start, end)) =
                        normalized_read_slice_bounds(read_start, read_end, read_len, strand)
                    {
                        assert!(end > start, "strand={strand}, start={start}, end={end}");
                        assert!(end <= read_len, "strand={strand}, end={end}, len={read_len}");
                    }
                }
            }
        }
    }
}
