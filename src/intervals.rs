use crate::extract;
use crate::methyl;
use crate::util;
use anyhow::{Context, Result as AnyhowResult};
use bio::io::fastq;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, warn};
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, Record};
use rust_htslib::faidx::Reader;
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
pub fn extract_haplotypes_coordinates_from_bam(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    primary_only: bool,
) -> AnyhowResult<(
    Vec<fastq::Record>,
    HashMap<String, (u64, u64)>,
    HashMap<String, String>,
    HashMap<String, Record>,
)> {
    let rg_sm_map = util::get_rg_to_sm_mapping(bam);
    let mut bmap: HashMap<String, (String, Vec<u8>)> = HashMap::new();
    let mut read_coordinates: HashMap<String, (u64, u64)> = HashMap::new(); // Track read coordinates
    let mut read_sequence_dict: HashMap<String, String> = HashMap::new();
    let mut bam_records: HashMap<String, Record> = HashMap::new();
    let mut eligible_reads: HashMap<String, (String, usize, usize)> = HashMap::new();
    // Create progress bar
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    pb.set_message("Scanning alignments...");

    let _ = bam.fetch((chr.as_bytes(), start, end));
    let mut seen_reads = HashSet::new();
    for record_result in bam.records() {
        let record = record_result?;
        let qname = String::from_utf8_lossy(record.qname()).into_owned();
        if seen_reads.contains(&qname) {
            continue;
        }
        seen_reads.insert(qname.clone());

        let is_secondary = record.is_secondary();
        let is_supplementary = record.is_supplementary();
        if primary_only && (is_secondary || is_supplementary) {
            continue;
        }

        let Some((query_start, query_end)) =
            extract::record_interval_query_bounds(&record, start, end)
        else {
            continue;
        };

        let sm = match util::get_sm_name_from_rg(&record, &rg_sm_map) {
            Ok(a) => a,
            Err(_) => String::from("unknown"),
        };
        let locus_key = format!("{}:{}-{}", chr, start, end);
        let seq_name = format!("{qname}|{locus_key}|{sm}");
        eligible_reads.insert(qname.clone(), (seq_name, query_start, query_end));
        read_sequence_dict.insert(
            qname.clone(),
            String::from_utf8_lossy(&record.seq().as_bytes()).into_owned(),
        );
        bam_records.insert(qname, record.clone());
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
                // skip the alignment from the same read
                if readnames.contains(&qname) {
                    continue;
                }
                readnames.insert(qname.clone());
                let Some((seq_name, _query_start, _query_end)) = eligible_reads.get(&qname) else {
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

    for (read_name, (_seq_name, query_start, query_end)) in eligible_reads.iter() {
        read_coordinates.insert(read_name.clone(), (*query_start as u64, *query_end as u64));
    }

    let records: Vec<fastq::Record> = eligible_reads
        .iter()
        .map(|(_read_name, (seq_name, _query_start, _query_end))| {
            let (seq, qual) = bmap
                .get(seq_name)
                .cloned()
                .unwrap_or_else(|| (String::new(), Vec::new()));
            fastq::Record::with_attrs(seq_name.as_str(), None, seq.as_bytes(), qual.as_slice())
        })
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
    let seq = fasta
        .fetch_seq_string(
            chr,
            usize::try_from(*start)?,
            usize::try_from(*stop - 1).unwrap(),
        )
        .unwrap();
    if !seq.is_empty() {
        let records = vec![fastq::Record::with_attrs(
            id.as_str(),
            None,
            seq.as_bytes(),
            vec![30; seq.len()].as_slice(),
        )];

        return Ok(records);
    }

    Err(anyhow::anyhow!("No sequence found for locus: {}", id))
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
    write_output: bool,
) -> AnyhowResult<(
    HashMap<String, (String, HashMap<String, HashMap<usize, f32>>, f64)>,
    HashMap<String, (u64, u64)>,
    HashMap<String, Vec<u8>>,
    HashMap<String, String>,
)> {
    let (
        read_coordinates,
        read_quality_dict,
        read_sequence_dict,
        read_strand_dict,
        bam_records_dict,
    ) = extract::extract_haplotypes_coordinates_from_bam(
        bam,
        chromosome,
        start as u64,
        end as u64,
        sampleid,
        primary_only,
    )
    .unwrap();
    let unique_local_sequences = read_sequence_dict.values().collect::<HashSet<_>>().len();
    debug!(
        "Window {}:{}-{} extracted reads: {}, unique local sequences: {}",
        chromosome,
        start,
        end,
        read_sequence_dict.len(),
        unique_local_sequences
    );
    let read_methyl_dict = methyl::start(bam_records_dict, &read_coordinates);
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
