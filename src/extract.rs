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

/// Unique key for one BAM alignment record (qname + 0-based reference start).
pub fn alignment_bounds_key(qname: &str, align_pos: i64) -> String {
    format!("{qname}|{align_pos}")
}

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


/// Per-alignment query bounds within a reference interval, plus window-edge alignment keys.
pub struct IntervalQueryBounds {
    /// First inclusive qpos in `[start, end)` keyed by `alignment_bounds_key`.
    pub per_align_start: HashMap<String, usize>,
    /// Last exclusive qpos in `[start, end)` keyed by `alignment_bounds_key`.
    pub per_align_end: HashMap<String, usize>,
    /// Alignment keys seen in the window with no valid qpos (fully deleted span).
    pub deletion_spanning: HashSet<String>,
    /// Read → alignment key active at the window start column.
    pub window_start_align: HashMap<String, String>,
    /// Read → alignment key active at the window end column.
    pub window_end_align: HashMap<String, String>,
}

/// Build the in-window read sequence, concatenating across supplementary alignments when
/// the window start and end fall on different alignment records for the same read.
pub fn sequence_from_interval_bounds(
    records: &[BamRecord],
    qname: &str,
    bounds: &IntervalQueryBounds,
) -> Option<String> {
    let start_key = bounds.window_start_align.get(qname)?;
    let end_key = bounds.window_end_align.get(qname)?;

    if start_key == end_key {
        let record = records
            .iter()
            .find(|r| alignment_bounds_key(qname, r.pos()) == *start_key)?;
        let start_pos = bounds.per_align_start.get(start_key)?;
        let end_pos = bounds.per_align_end.get(end_key)?;
        let (pos_start, pos_end) =
            normalized_read_slice_bounds(*start_pos, *end_pos, record.seq().len(), "")?;
        return Some(String::from_utf8_lossy(&record.seq().as_bytes()[pos_start..pos_end]).to_string());
    }

    let start_ref_pos: i64 = start_key
        .rsplit('|')
        .next()?
        .parse()
        .ok()?;
    let end_ref_pos: i64 = end_key.rsplit('|').next()?.parse().ok()?;

    let mut segments: Vec<(i64, String)> = Vec::new();
    for record in records {
        if String::from_utf8_lossy(record.qname()) != qname {
            continue;
        }
        let align_key = alignment_bounds_key(qname, record.pos());
        let ref_pos = record.pos();
        if ref_pos < start_ref_pos || ref_pos > end_ref_pos {
            continue;
        }
        let Some(start_pos) = bounds.per_align_start.get(&align_key) else {
            continue;
        };
        let Some(end_pos) = bounds.per_align_end.get(&align_key) else {
            continue;
        };
        let Some((pos_start, pos_end)) =
            normalized_read_slice_bounds(*start_pos, *end_pos, record.seq().len(), "")
        else {
            continue;
        };
        let seq = String::from_utf8_lossy(&record.seq().as_bytes()[pos_start..pos_end]).to_string();
        segments.push((ref_pos, seq));
    }
    segments.sort_by_key(|(pos, _)| *pos);
    if segments.is_empty() {
        return None;
    }
    Some(segments.into_iter().map(|(_, seq)| seq).collect())
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
) -> AnyhowResult<IntervalQueryBounds> {
    let mut per_align_start: HashMap<String, usize> = HashMap::new();
    let mut per_align_end: HashMap<String, usize> = HashMap::new();
    let mut seen_in_window: HashSet<String> = HashSet::new();
    let mut window_start_align: HashMap<String, String> = HashMap::new();
    let mut window_end_align: HashMap<String, String> = HashMap::new();

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
            let align_key = alignment_bounds_key(&qname, record.pos());
            if seen.contains(&align_key) {
                continue;
            }
            seen.insert(align_key.clone());
            seen_in_window.insert(align_key.clone());
            if let Some(qpos) = alignment.qpos() {
                per_align_start.entry(align_key.clone()).or_insert(qpos);
                let excl = match alignment.indel() {
                    bam::pileup::Indel::Ins(len) => qpos + len as usize + 1,
                    bam::pileup::Indel::Del(_) | bam::pileup::Indel::None => qpos + 1,
                };
                per_align_end.insert(align_key.clone(), excl);
                if ref_pos == start {
                    window_start_align
                        .entry(qname.clone())
                        .or_insert(align_key.clone());
                }
                window_end_align.insert(qname, align_key);
            }
        }
    }

    let deletion_spanning: HashSet<String> = seen_in_window
        .into_iter()
        .filter(|key| !per_align_start.contains_key(key))
        .collect();

    Ok(IntervalQueryBounds {
        per_align_start,
        per_align_end,
        deletion_spanning,
        window_start_align,
        window_end_align,
    })
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
    let bounds = record_interval_query_bounds_range(bam, start, end)?;

    let mut filtered_sequence_dict = HashMap::new();
    let _ = bam.fetch((chr.as_bytes(), start, end));
    let mut records_by_qname: HashMap<String, Vec<BamRecord>> = HashMap::new();
    for record_result in bam.records() {
        let record = record_result?;
        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            continue;
        }
        let read_name = String::from_utf8_lossy(record.qname()).into_owned();
        records_by_qname
            .entry(read_name)
            .or_default()
            .push(record);
    }

    for (read_name, mut records) in records_by_qname {
        records.sort_by_key(|r| r.pos());
        let start_key = match bounds.window_start_align.get(&read_name) {
            Some(k) => k.clone(),
            None => continue,
        };
        if bounds.deletion_spanning.contains(&start_key) {
            let record_id = format!("{read_name}|{chr}:{start}-{end}|{sampleid}");
            filtered_sequence_dict.insert(record_id, String::new());
            continue;
        }

        let Some(read_seq_final) = sequence_from_interval_bounds(&records, &read_name, &bounds) else {
            continue;
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
    let bounds = record_interval_query_bounds_range(bam, start, end)?;

    let mut read_coordinates_formatted = HashMap::new();
    let mut read_sequence_dict_formatted = HashMap::new();
    let mut read_quality_dict_formatted = HashMap::new();
    let mut read_strand_dict_formatted = HashMap::new();
    let mut bam_records_dict_formatted = HashMap::new();
    let _ = bam.fetch((chr.as_bytes(), start, end));
    for record_result in bam.records() {
        let record = record_result?;
        let read_name = String::from_utf8_lossy(record.qname()).into_owned();
        let align_key = alignment_bounds_key(&read_name, record.pos());

        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            continue;
        }

        if bounds.deletion_spanning.contains(&align_key) {
            let record_id = format!(
                "{}|{}",
                intervals::generate_read_name(&read_name, chr, start, end, sampleid),
                record.pos()
            );
            read_coordinates_formatted.insert(record_id.clone(), (0u64, 0u64));
            read_sequence_dict_formatted.insert(record_id.clone(), Vec::new());
            read_quality_dict_formatted.insert(record_id.clone(), Vec::new());
            read_strand_dict_formatted
                .insert(record_id.clone(), record.strand().to_string());
            bam_records_dict_formatted.insert(record_id, record.clone());
            continue;
        }

        let Some(start_pos) = bounds.per_align_start.get(&align_key) else {
            continue;
        };
        let Some(end_pos) = bounds.per_align_end.get(&align_key) else {
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

        let read_seq_final = String::from_utf8_lossy(&read_seq[pos_start..pos_end]).to_string();

        let record_id = format!(
            "{}|{}",
            intervals::generate_read_name(&read_name, chr, start, end, sampleid),
            record.pos()
        );
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
