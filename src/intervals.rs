use crate::extract;
use crate::methyl;
use crate::util;
use anyhow::{Context, Result as AnyhowResult};
use bio::io::fastq;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, warn};
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, Record};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

pub fn generate_read_name(
    read_name: &String,
    chr: &str,
    start: u64,
    end: u64,
    sampleid: &String,
) -> String {
    format!("{read_name}|{chr}:{start}-{end}|{sampleid}")
}

// extract haplotypes from bam file
pub fn extract_haplotypes_coordinates_from_bam_pileup(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    primary_only: bool,
) -> AnyhowResult<(
    Vec<fastq::Record>,
    HashMap<String, (u64, u64)>,
    HashMap<String, Vec<u8>>,
    HashMap<String, HashMap<usize, f32>>,
    HashMap<String, String>,
)> {
    let rg_sm_map = util::get_rg_to_sm_mapping(bam);
    let mut bmap: HashMap<String, (String, Vec<u8>)> = HashMap::new();
    let mut read_coordinates: HashMap<String, (u64, u64)> = HashMap::new(); // Track read coordinates
    let mut read_sequence_dict: HashMap<String, Vec<u8>> = HashMap::new();
    let mut read_methyl_dict: HashMap<String, HashMap<usize, f32>> = HashMap::new();
    let mut read_strand_dict: HashMap<String, String> = HashMap::new();
    let mut eligible_reads: HashMap<String, (String, usize, usize)> = HashMap::new();
    // Create progress bar
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    // pb.set_message("Scanning alignments...");

    let _ = bam.fetch((chr.as_bytes(), start, end));
    let mut collected_records: Vec<Record> = Vec::new();
    let mut seen_reads = HashSet::new();
    for record_result in bam.records() {
        let record = record_result?;
        let qname = String::from_utf8_lossy(record.qname()).into_owned();
        let align_key = extract::alignment_bounds_key(&qname, record.pos());
        if seen_reads.contains(&align_key) {
            continue;
        }
        seen_reads.insert(align_key.clone());
        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            continue;
        }
        collected_records.push(record);
    }

    let _ = bam.fetch((chr.as_bytes(), start, end));
    let bounds = extract::record_interval_query_bounds_range(bam, start, end)?;

    for record in collected_records {
        let qname = String::from_utf8_lossy(record.qname()).into_owned();
        let align_key = extract::alignment_bounds_key(&qname, record.pos());
        let sm = match util::get_sm_name_from_rg(&record, &rg_sm_map) {
            Ok(a) => a,
            Err(_) => String::from("unknown"),
        };
        let locus_key = format!("{}:{}-{}", chr, start, end);
        let seq_name = format!("{qname}|{locus_key}|{sm}");

        if bounds.deletion_spanning.contains(&align_key)
            && !bounds.per_align_start.contains_key(&align_key)
        {
            eligible_reads.insert(align_key.clone(), (seq_name.clone(), 0, 0));
            bmap.insert(seq_name, (String::new(), Vec::new()));
            read_strand_dict.insert(align_key.clone(), record.strand().to_string());
            read_methyl_dict.insert(align_key, HashMap::new());
            continue;
        }

        let window_start_key = bounds.window_start_align.get(&qname);
        let window_end_key = bounds.window_end_align.get(&qname);
        let query_start = window_start_key
            .and_then(|k| bounds.per_align_start.get(k))
            .copied();
        let query_end = window_end_key
            .and_then(|k| bounds.per_align_end.get(k))
            .copied();
        let (Some(query_start), Some(query_end)) = (query_start, query_end) else {
            continue;
        };

        eligible_reads.insert(align_key.clone(), (seq_name, query_start, query_end));
        read_strand_dict.insert(align_key.clone(), record.strand().to_string());
        read_methyl_dict.insert(
            align_key,
            methyl::get_methylation_read(&record, query_start, query_end, 'm'),
        );
    }

    let _ = bam.fetch((chr.as_bytes(), start, end));
    pb.set_message("Processing pileup...");
    for p in bam.pileup() {
        let pileup = p?;
        pb.set_message(format!("Processing position {}", pileup.pos()));
        let mut readnames = HashSet::new();
        if start <= (pileup.pos() as u64) && (pileup.pos() as u64) < end {
            for alignment in pileup.alignments() {
                let record = alignment.record();
                let qname = String::from_utf8_lossy(record.qname()).into_owned();
                let align_key = extract::alignment_bounds_key(&qname, record.pos());
                if readnames.contains(&align_key) {
                    continue;
                }
                readnames.insert(align_key.clone());
                let Some((seq_name, _query_start, _query_end)) = eligible_reads.get(&align_key) else {
                    continue;
                };

                if !bmap.contains_key(seq_name) {
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
                                bmap.get_mut(seq_name).unwrap().0.push(a as char);
                                bmap.get_mut(seq_name).unwrap().1.push(valid_q + 33);
                            }
                        }
                    }
                    bam::pileup::Indel::Del(_) => {
                        // For deletions, add the first base of the deletion
                        if let Some(qpos) = alignment.qpos() {
                            let a = record.seq()[qpos];
                            let q = record.qual()[qpos];
                            let valid_q = q.min(40);

                            bmap.get_mut(seq_name).unwrap().0.push(a as char);
                            bmap.get_mut(seq_name).unwrap().1.push(valid_q + 33);
                        }
                    }
                    bam::pileup::Indel::None => {
                        // For matches/mismatches, add the base

                        if let Some(qpos) = alignment.qpos() {
                            let a = record.seq()[qpos];
                            let q = record.qual()[qpos];
                            let valid_q = q.min(40);

                            bmap.get_mut(seq_name).unwrap().0.push(a as char);
                            bmap.get_mut(seq_name).unwrap().1.push(valid_q + 33);
                        }
                    }
                }
            }
        }
    }

    pb.finish_with_message("Pileup processing completed");

    read_sequence_dict.clear();
    for (read_name, (seq_name, query_start, query_end)) in eligible_reads.iter() {
        read_coordinates.insert(read_name.clone(), (*query_start as u64, *query_end as u64));
        // Include deletion-spanning reads: their bmap entry has an empty sequence, which is
        // intentional — they map to a large deletion covering the entire window.
        if let Some((seq, _qual)) = bmap.get(seq_name) {
            read_sequence_dict.insert(read_name.clone(), seq.as_bytes().to_vec());
        }
    }

    let mut seen_seq_names = HashSet::new();
    let records: Vec<fastq::Record> = eligible_reads
        .iter()
        .filter_map(|(_read_name, (seq_name, _query_start, _query_end))| {
            if !seen_seq_names.insert(seq_name.clone()) {
                return None;
            }
            let (seq, qual) = bmap.get(seq_name)?;
            Some(fastq::Record::with_attrs(
                seq_name.as_str(),
                None,
                seq.as_bytes(),
                qual.as_slice(),
            ))
        })
        .collect();

    debug!(
        "Extracted {} reads from region ({} with pileup sequences)",
        eligible_reads.len(),
        read_sequence_dict.len()
    );

    Ok((records, read_coordinates, read_sequence_dict, read_methyl_dict, read_strand_dict))
}

pub fn process_fasta_file(
    reference: &Vec<fastq::Record>,
    chromosome: &str,
    start: usize,
    end: usize,
    sampleid: &String,
) -> Vec<fastq::Record> {
    let mut reference_seqs = Vec::new();
    for record in reference.iter() {
        let seq_id = record.id().to_string();
        let sequence = String::from_utf8_lossy(record.seq()).to_string();
        // Check if this sequence matches the target chromosome
        if seq_id == chromosome {
            debug!(
                "Extracting region {}:{}-{} from chromosome {}",
                chromosome, start, end, seq_id
            );
            // Extract the specified region
            if start < sequence.len() && end <= sequence.len() {
                let region_seq = sequence[start..end].to_string();
                let record_id = format!("{}:{}-{}|reference|{}", chromosome, start, end, sampleid);

                let fastq_record = fastq::Record::with_attrs(
                    &record_id,
                    None,
                    region_seq.as_bytes(),
                    vec![30; region_seq.len()].as_slice(), // Default quality score
                );

                reference_seqs.push(fastq_record);
                debug!("Extracted reference sequence: {} bp", region_seq.len());
            } else {
                warn!(
                    "Region coordinates out of bounds for sequence length {}",
                    sequence.len()
                );
            }
        }
    }

    if reference_seqs.is_empty() {
        warn!(
            "No matching chromosome '{}' found in FASTA file",
            chromosome
        );
    }

    debug!("Extracted {} reference sequences", reference_seqs.len());
    reference_seqs
}

pub fn write_fasta_output(
    final_hap: HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)>,
    output_path: &PathBuf,
) -> AnyhowResult<()> {
    let mut file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;

    let mut index: i32 = 1;
    let chars_per_line = 60;

    for (sequence, (cigar, read_names, allele_frequency)) in final_hap {
        // Use the first read name as the sequence ID, or create a hash-based ID
        let sequence_id = if !read_names.is_empty() {
            serde_json::to_string(&read_names).expect("Failed to serialize to JSON")
        } else {
            "unknown".to_string()
        };
        writeln!(
            file,
            ">Haplotype {}\tLength{}\tCigar {}\tSupportReads {}\tAlleleFrequency {:.2}\t{}",
            index,
            sequence.len(),
            cigar,
            read_names.len(),
            allele_frequency,
            sequence_id
        )?;
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

pub fn collapse_haplotypes(
    read_seq_dict: &HashMap<String, Vec<u8>>,
    read_methyl_dict: &HashMap<String, HashMap<usize, f32>>,
    reference: &fastq::Record,
    min_reads: usize,
    frequency_min: f64,
) -> AnyhowResult<HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)>> {
    // Read the reads records (name and sequence) into a vector.
    let mut read_dictionary: HashMap<String, HashMap<String, HashMap<usize, f32>>> = HashMap::new();
    for (read_name, read_seq) in read_seq_dict.iter() {
        let read_name_origin = read_name.to_string();
        let contig = String::from_utf8_lossy(read_seq).to_string();
        if !read_methyl_dict.contains_key(&read_name_origin) {
            read_dictionary
                .entry(contig)
                .or_default()
                .insert(read_name_origin, HashMap::new());
        } else {
            let methyl_dict = read_methyl_dict.get(&read_name_origin).unwrap();
            read_dictionary
                .entry(contig)
                .or_default()
                .insert(read_name_origin, methyl_dict.clone());
        }
    }

    let reference_seq = String::from_utf8_lossy(reference.seq()).to_string();
    let mut final_hap: HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)> =
        HashMap::new();

    for (hap, vec) in read_dictionary.iter() {
        let allele_frequency = vec.len() as f64 / read_seq_dict.len() as f64;
        if (vec.len() >= min_reads) && (allele_frequency >= frequency_min) {
            let cigar = util::gap_open_aligner(&reference_seq, hap);
            final_hap.insert(hap.clone(), (cigar, vec.clone(), allele_frequency));
        }
    }
    debug!(
        "Extracted {} haplotypes from region, total reads: {:?}",
        final_hap.len(),
        read_dictionary
            .values()
            .map(|x| x.len())
            .collect::<Vec<_>>()
    );
    Ok(final_hap)
}

pub fn start(
    bam: &mut IndexedReader,
    reference_fa: &Vec<fastq::Record>,
    chromosome: &str,
    start: usize,
    end: usize,
    sampleid: &String,
    min_reads: usize,
    frequency_min: f64,
    primary_only: bool,
    pileup:bool,
    write_output: bool,
) -> AnyhowResult<(
    HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)>,
    HashMap<String, (u64, u64)>,
    HashMap<String, Vec<u8>>,
    HashMap<String, String>,
)> {
    // if pileup is true, use extract_haplotypes_coordinates_from_bam_pileup
    // if pileup is false, use  extract::extract_haplotypes_coordinates_from_bam
    let (read_coordinates, read_quality_dict, read_sequence_dict, read_strand_dict, read_methyl_dict) = if !pileup {
        extract::extract_haplotypes_coordinates_from_bam(
            bam,
            chromosome,
            start as u64,
            end as u64,
            sampleid,
            primary_only,
        )
        .unwrap()
    } else {
        let (reads, read_coordinates, read_sequence_dict, read_methyl_dict, read_strand_dict) =
            extract_haplotypes_coordinates_from_bam_pileup(
                bam,
                chromosome,
                start as u64,
                end as u64,
                primary_only,
            )
            .unwrap();
        let _ = reads;
        (
            read_coordinates,
            HashMap::<String, Vec<u8>>::new(),
            read_sequence_dict,
            read_strand_dict,
            read_methyl_dict,
        )
    };
    // let (reads, read_coordinates, read_sequence_dict, bam_records_dict, read_strand_dict) =
    // extract_haplotypes_coordinates_from_bam_pileup(
    //     bam,
    //     chromosome,
    //     start as u64,
    //     end as u64,
    //     primary_only,
    // )
    // .unwrap();
    let unique_local_sequences = read_sequence_dict.values().collect::<HashSet<_>>().len();
    debug!(
        "Window {}:{}-{} extracted reads: {}, unique local sequences: {}",
        chromosome,
        start,
        end,
        read_sequence_dict.len(),
        unique_local_sequences
    );
    let _ = read_quality_dict;
    debug!("Extracted {} reads from region", read_coordinates.len());

    // let reads_list = process_bam_file_by_coordinates(bam, chromosome, start, end, primary_only, sampleid);
    let reference = process_fasta_file(reference_fa, chromosome, start, end, sampleid);
    let reference = reference.first().unwrap().clone();
    let final_hap = collapse_haplotypes(
        &read_sequence_dict,
        &read_methyl_dict,
        &reference,
        min_reads,
        frequency_min,
    )?;
    debug!(
        "Window {}:{}-{} retained haplotypes after filtering: {} (min_reads={}, frequency_min={})",
        chromosome,
        start,
        end,
        final_hap.len(),
        min_reads,
        frequency_min
    );
    debug!("Haplotype reconstruction completed");

    if write_output {
        let final_hap_output = PathBuf::from(format!(
            "{}/{}_{}_{}_{}_haplograph.fasta",
            ".", sampleid, chromosome, start, end
        ));
        write_fasta_output(final_hap.clone(), &final_hap_output)?;
    }

    Ok((
        final_hap,
        read_coordinates,
        read_sequence_dict,
        read_strand_dict,
    ))
}

#[cfg(test)]
mod pileup_extract_tests {
    use super::{collapse_haplotypes, extract_haplotypes_coordinates_from_bam_pileup, process_fasta_file};
    use bio::io::fastq;
    use rust_htslib::bam::IndexedReader;
    use std::collections::HashSet;

    fn chr17_reference(len: usize) -> Vec<fastq::Record> {
        vec![fastq::Record::with_attrs(
            "chr17",
            None,
            vec![b'N'; len].as_slice(),
            &[],
        )]
    }

    #[test]
    fn pileup_extraction_produces_collapsible_haplotypes() {
        let bam_path = "HG00097_tp53_test.bam";
        if !std::path::Path::new(bam_path).exists() {
            return;
        }

        let chr = "chr17";
        let start = 7668752_u64;
        let end = 7668828_u64;
        let sampleid = "HG00097".to_string();

        let mut bam = IndexedReader::from_path(bam_path).unwrap();
        let (records, _coords, read_seq_dict, read_methyl, _strand) =
            extract_haplotypes_coordinates_from_bam_pileup(&mut bam, chr, start, end, false)
                .unwrap();

        let unique_pileup: HashSet<_> = records
            .iter()
            .map(|r| String::from_utf8_lossy(r.seq()).into_owned())
            .collect();
        let unique_for_collapse: HashSet<_> = read_seq_dict
            .values()
            .map(|s| String::from_utf8_lossy(s).into_owned())
            .collect();

        assert_eq!(records.len(), read_seq_dict.len());
        assert_eq!(unique_pileup.len(), unique_for_collapse.len());
        assert!(unique_for_collapse.len() <= 10, "got {}", unique_for_collapse.len());

        let reference =
            process_fasta_file(&chr17_reference(end as usize + 1), chr, start as usize, end as usize, &sampleid)
                .first()
                .unwrap()
                .clone();
        let final_hap = collapse_haplotypes(&read_seq_dict, &read_methyl, &reference, 2, 0.0).unwrap();

        assert!(
            !final_hap.is_empty(),
            "expected haplotype nodes after collapse, unique_seqs={}",
            unique_for_collapse.len()
        );
    }
}
