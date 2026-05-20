use crate::intervals;
use crate::util;
use anyhow::{Context, Result as AnyhowResult};
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use log::info;
use rust_htslib::bam::{record::Cigar, IndexedReader, Read as BamRead, Record as BamRecord};
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

pub(crate) fn query_boundary_offset_from_cigar(
    cigar: &[Cigar],
    ref_start: u64,
    target_ref: u64,
) -> Option<usize> {
    let mut ref_pos = ref_start;
    let mut query_pos = 0_usize;

    for op in cigar {
        match *op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = len as u64;
                if target_ref >= ref_pos && target_ref <= ref_pos + len {
                    return Some(query_pos + (target_ref - ref_pos) as usize);
                }
                ref_pos += len;
                query_pos += len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                let len = len as u64;
                if target_ref >= ref_pos && target_ref <= ref_pos + len {
                    return Some(query_pos);
                }
                ref_pos += len;
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                query_pos += len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    if target_ref == ref_pos {
        Some(query_pos)
    } else {
        None
    }
}

pub(crate) fn record_interval_query_bounds(
    record: &BamRecord,
    start: u64,
    end: u64,
) -> Option<(usize, usize)> {
    if start > end {
        return None;
    }

    let ref_start = record.pos() as u64;
    let ref_end = record.cigar().end_pos() as u64;
    if start < ref_start || end > ref_end {
        return None;
    }

    let cigar = record.cigar().iter().copied().collect::<Vec<_>>();
    let query_start = query_boundary_offset_from_cigar(&cigar, ref_start, start)?;
    let query_end = query_boundary_offset_from_cigar(&cigar, ref_start, end)?;
    Some((query_start, query_end))
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
    if pos_end < pos_start || pos_end > read_len {
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
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    // pb.set_message("Processing alignments...");

    let _ = bam.fetch((chr.as_bytes(), start, end));

    let mut filtered_sequence_dict = HashMap::new();
    let mut seen_reads = HashSet::new();
    for record_result in bam.records() {
        let record = record_result?;
        let read_name = String::from_utf8_lossy(record.qname()).into_owned();
        if seen_reads.contains(&read_name) {
            continue;
        }
        seen_reads.insert(read_name.clone());

        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            continue;
        }

        let Some((read_start, read_end)) = record_interval_query_bounds(&record, start, end) else {
            continue;
        };

        let read_seq = record.seq().as_bytes();
        let read_strand = record.strand().to_string();
        let read_len = read_seq.len();
        let Some((pos_start, pos_end)) =
            normalized_read_slice_bounds(read_start, read_end, read_len, &read_strand)
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
        filtered_sequence_dict.insert(record_id, read_seq_final);
        if filtered_sequence_dict.len() % 100 == 0 {
            pb.set_message(format!("Collected {} reads", filtered_sequence_dict.len()));
        }
    }

    pb.finish_with_message("Alignment processing completed");
    info!(
        "Extracted {} reads from region",
        filtered_sequence_dict.len()
    );

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
) -> AnyhowResult<(
    HashMap<String, (u64, u64)>,
    HashMap<String, Vec<u8>>,
    HashMap<String, Vec<u8>>,
    HashMap<String, String>,
    HashMap<String, BamRecord>,
)> {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    // pb.set_message("Processing alignments...");

    let _ = bam.fetch((chr.as_bytes(), start, end));
    let mut read_coordinates_formatted = HashMap::new();
    let mut read_sequence_dict_formatted = HashMap::new();
    let mut read_quality_dict_formatted = HashMap::new();
    let mut read_strand_dict_formatted = HashMap::new();
    let mut bam_records_dict_formatted = HashMap::new();
    let mut seen_reads = HashSet::new();
    for record_result in bam.records() {
        let record = record_result?;
        let read_name = String::from_utf8_lossy(record.qname()).into_owned();
        if seen_reads.contains(&read_name) {
            continue;
        }
        seen_reads.insert(read_name.clone());

        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            continue;
        }

        let Some((read_start, read_end)) = record_interval_query_bounds(&record, start, end) else {
            continue;
        };

        let read_seq = record.seq().as_bytes();
        let read_strand = record.strand().to_string();
        let read_len = read_seq.len();
        let Some((pos_start, pos_end)) =
            normalized_read_slice_bounds(read_start, read_end, read_len, &read_strand)
        else {
            continue;
        };

        let read_seq_sub = String::from_utf8_lossy(&read_seq[pos_start..pos_end]).to_string();
        let read_seq_final = if read_strand == "+" {
            read_seq_sub
        } else {
            util::reverse_complement(&read_seq_sub)
        };

        let record_id = intervals::generate_read_name(&read_name, chr, start, end, sampleid);
        read_coordinates_formatted.insert(record_id.clone(), (read_start as u64, read_end as u64));
        read_sequence_dict_formatted.insert(record_id.clone(), read_seq_final.into_bytes());
        read_quality_dict_formatted.insert(record_id.clone(), record.qual().to_vec());
        read_strand_dict_formatted.insert(record_id.clone(), read_strand);
        bam_records_dict_formatted.insert(record_id, record.clone());

        if read_coordinates_formatted.len() % 100 == 0 {
            pb.set_message(format!(
                "Collected {} reads",
                read_coordinates_formatted.len()
            ));
        }
    }

    pb.finish_with_message("Alignment processing completed");
    info!(
        "Extracted {} reads from region",
        read_coordinates_formatted.len()
    );

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
    pileup: bool,
) -> AnyhowResult<()> {
    if pileup {
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
    } else {
        let reads = extract_haplotypes_from_bam(
            bam,
            &chromosome,
            start as u64,
            end as u64,
            &sampleid,
            primary_only,
        )
        .unwrap();

        let output_p = PathBuf::from(format!("{}.fasta", output_path));
        if !reads.is_empty() {
            let _ = util::write_fasta(&reads, &output_p);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{normalized_read_slice_bounds, query_boundary_offset_from_cigar};
    use rust_htslib::bam::record::Cigar;

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
    fn reverse_strand_equal_coordinates_preserve_empty_deleted_window() {
        let bounds = normalized_read_slice_bounds(42, 42, 100, "-").unwrap();
        assert_eq!(bounds, (58, 58));
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
                        assert!(end >= start, "strand={strand}, start={start}, end={end}");
                        assert!(
                            end <= read_len,
                            "strand={strand}, end={end}, len={read_len}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn cigar_boundary_inside_deletion_returns_same_query_offset() {
        let cigar = vec![Cigar::Match(50), Cigar::Del(200), Cigar::Match(50)];
        let start = query_boundary_offset_from_cigar(&cigar, 100, 200).unwrap();
        let end = query_boundary_offset_from_cigar(&cigar, 100, 300).unwrap();
        assert_eq!((start, end), (50, 50));
    }

    #[test]
    fn cigar_boundary_after_deletion_resumes_query_progress() {
        let cigar = vec![Cigar::Match(50), Cigar::Del(200), Cigar::Match(50)];
        let start = query_boundary_offset_from_cigar(&cigar, 100, 300).unwrap();
        let end = query_boundary_offset_from_cigar(&cigar, 100, 380).unwrap();
        assert_eq!((start, end), (50, 80));
    }
}
