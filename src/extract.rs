use crate::intervals;
use crate::util;
use anyhow::{Context, Result as AnyhowResult};
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use log::info;
use rust_htslib::bam::{IndexedReader, Read as BamRead, Record as BamRecord};
use rust_htslib::bam;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;


pub fn normalized_read_slice_bounds(
    read_start: usize,
    read_end: usize,
    read_len: usize,
    _read_strand: &str,
) -> Option<(usize, usize)> {
    // BAM always stores sequences in forward-strand orientation for both plus and minus
    // strand reads (FLAG 0x10 means the stored SEQ is already the reverse-complement of
    // the original read, i.e. it is in the reference left→right direction).  Therefore
    // qpos always increases with reference position for both strands, start_pos < end_pos
    // in all normal cases, and no strand-specific coordinate transformation is needed.
    let pos_start = read_start.min(read_end);
    let pos_end = read_start.max(read_end);
    if pos_end == pos_start || pos_end > read_len {
        return None;
    }
    Some((pos_start, pos_end))
}


/// Single-pass pileup scan over `[start, end)`.
///
/// Returns `(start_dict, end_dict, deletion_spanning)` where:
///   - `start_dict`: each read's **first** valid inclusive qpos in the window.
///     Reads whose base at exactly `start` is deleted are still included using the
///     first non-deleted position they have within `[start, end)`.
///   - `end_dict`: each read's **last** valid exclusive qpos in the window.
///     Reads whose base at exactly `end-1` is deleted are still included using the
///     last non-deleted position they have within `[start, end)`.
///   - `deletion_spanning`: reads that appear in the pileup within `[start, end)` but
///     have `qpos = None` for every column — i.e. the entire query window falls inside
///     a large deletion of that read.  Callers should emit these reads with an empty
///     sequence (`""`).
///
/// Exclusive-end semantics per indel type at the last valid column:
///   - `None`   → `qpos + 1`         (include the matched base)
///   - `Ins(N)` → `qpos + N + 1`     (include anchor base + inserted bases)
///   - `Del(_)` → `qpos + 1`         (deletion starts outside the window; include the base)
pub fn record_interval_query_bounds_range(
    bam: &mut IndexedReader,
    start: u64,
    end: u64,
) -> AnyhowResult<(HashMap<String, usize>, HashMap<String, usize>, HashSet<String>)> {
    let mut start_dict: HashMap<String, usize> = HashMap::new();
    let mut end_dict: HashMap<String, usize> = HashMap::new();
    // Track every read that appears in any pileup column within [start, end).
    // Reads here but absent from start_dict have qpos=None for the entire window
    // (i.e. the window lies completely inside a deletion of that read).
    let mut seen_in_window: HashSet<String> = HashSet::new();
    for p in bam.pileup() {
        let pileup = p?;
        let ref_pos = pileup.pos() as u64;
        if ref_pos < start {
            continue;
        }
        if ref_pos >= end {
            break;
        }
        let mut seen: HashSet<String> = HashSet::new();
        for alignment in pileup.alignments() {
            let record = alignment.record();
            let qname = String::from_utf8_lossy(record.qname()).into_owned();
            if seen.contains(&qname) {
                continue;
            }
            seen.insert(qname.clone());
            seen_in_window.insert(qname.clone());
            if let Some(qpos) = alignment.qpos() {
                // First valid position → start boundary (keep first, never overwrite)
                start_dict.entry(qname.clone()).or_insert(qpos);
                // Last valid exclusive position → end boundary (always overwrite to keep last)
                let excl = match alignment.indel() {
                    bam::pileup::Indel::Ins(len) => qpos + len as usize + 1,
                    bam::pileup::Indel::Del(_) | bam::pileup::Indel::None => qpos + 1,
                };
                end_dict.insert(qname, excl);
            }
        }
    }
    // Reads seen in the window but with no valid qpos in any column are entirely
    // covered by a deletion — return them separately so callers can emit "".
    let deletion_spanning: HashSet<String> = seen_in_window
        .into_iter()
        .filter(|name| !start_dict.contains_key(name))
        .collect();
    Ok((start_dict, end_dict, deletion_spanning))
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
    let (start_pos_dict, end_pos_dict, deletion_spanning) =
        record_interval_query_bounds_range(bam, start, end)?;

    let mut filtered_sequence_dict = HashMap::new();
    let mut seen_reads = HashSet::new();
    let _ = bam.fetch((chr.as_bytes(), start, end));
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

        // The entire window lies within a deletion of this read → emit with empty sequence.
        if deletion_spanning.contains(&read_name) {
            let record_id = format!("{read_name}|{chr}:{start}-{end}|{sampleid}");
            filtered_sequence_dict.insert(record_id, String::new());
            continue;
        }

        let Some(start_pos) = start_pos_dict.get(&read_name) else {
            continue;
        };
        let Some(end_pos) = end_pos_dict.get(&read_name) else {
            continue;
        };

        let read_seq = record.seq().as_bytes();
        let read_strand = record.strand().to_string();
        let read_len = read_seq.len();
        let Some((pos_start, pos_end)) =
            normalized_read_slice_bounds(*start_pos, *end_pos, read_len, &read_strand)
        else {
            continue;
        };

        // BAM stores sequences in forward-strand orientation for both strands; no RC needed.
        let read_seq_final = String::from_utf8_lossy(&read_seq[pos_start..pos_end]).to_string();

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
    let (start_pos_dict, end_pos_dict, deletion_spanning) =
        record_interval_query_bounds_range(bam, start, end)?;

    let mut read_coordinates_formatted = HashMap::new();
    let mut read_sequence_dict_formatted = HashMap::new();
    let mut read_quality_dict_formatted = HashMap::new();
    let mut read_strand_dict_formatted = HashMap::new();
    let mut bam_records_dict_formatted = HashMap::new();
    let mut seen_reads = HashSet::new();
    let _ = bam.fetch((chr.as_bytes(), start, end));
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

        // The entire window lies within a deletion of this read → emit with empty sequence.
        if deletion_spanning.contains(&read_name) {
            let record_id = intervals::generate_read_name(&read_name, chr, start, end, sampleid);
            read_coordinates_formatted.insert(record_id.clone(), (0u64, 0u64));
            read_sequence_dict_formatted.insert(record_id.clone(), Vec::new());
            read_quality_dict_formatted.insert(record_id.clone(), Vec::new());
            read_strand_dict_formatted
                .insert(record_id.clone(), record.strand().to_string());
            bam_records_dict_formatted.insert(record_id, record.clone());
            continue;
        }

        let Some(start_pos) = start_pos_dict.get(&read_name) else {
            continue;
        };
        let Some(end_pos) = end_pos_dict.get(&read_name) else {
            continue;
        };
        let read_seq = record.seq().as_bytes();
        let read_strand = record.strand().to_string();
        let read_len = read_seq.len();
        let Some((pos_start, pos_end)) =
            normalized_read_slice_bounds(*start_pos, *end_pos, read_len, &read_strand)
        else {
            continue;
        };

        // BAM stores sequences in forward-strand orientation for both strands; no RC needed.
        let read_seq_final = String::from_utf8_lossy(&read_seq[pos_start..pos_end]).to_string();

        let record_id = intervals::generate_read_name(&read_name, chr, start, end, sampleid);
        read_coordinates_formatted.insert(record_id.clone(), (*start_pos as u64, *end_pos as u64));
        read_sequence_dict_formatted.insert(record_id.clone(), read_seq_final.clone().into_bytes());
        read_quality_dict_formatted.insert(record_id.clone(), record.qual().to_vec());
        read_strand_dict_formatted.insert(record_id.clone(), read_strand);
        bam_records_dict_formatted.insert(record_id.clone(), record.clone());

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
        let (reads, read_coordinates, read_sequence_dictionary, bam_records, read_strand_dictionary) =
            intervals::extract_haplotypes_coordinates_from_bam_pileup(
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
mod extract_tests {
    use super::extract_haplotypes_from_bam;
    use rust_htslib::bam::IndexedReader;

    #[test]
    fn extract_haplotypes_from_bam_tp53_interval() {
        let bam_path = "HG00097_tp53_test.bam";
        if !std::path::Path::new(bam_path).exists() {
            return;
        }
        let mut bam = IndexedReader::from_path(bam_path).unwrap();
        let sampleid = "HG00097".to_string();
        let reads = extract_haplotypes_from_bam(
            &mut bam,
            "chr17",
            7668752,
            7668828,
            &sampleid,
            false,
        )
        .unwrap();
        assert_eq!(reads.len(), 34);
    }
}
