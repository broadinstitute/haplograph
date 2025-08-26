use log::{info, warn};
use anyhow::{Result, Context};
use rust_htslib::bam::IndexedReader;
use std::path::{PathBuf};
use url::Url;
use bio::io::fastq;
use std::collections::{HashMap, BTreeMap};
use std::io::Write;
use std::fs::File;
use crate::Cli;
use crate::util;
use std::error::Error;

pub fn open_bam_file(cli: &Cli) -> IndexedReader {
        // Open BAM file
        if cli.alignment_bam.to_string().starts_with("gs://") && std::env::var("GCS_OAUTH_TOKEN").is_err() {
            util::gcs_authorize_data_access();
        }
        let mut bam = if cli.alignment_bam.starts_with("gs://") {
            let url = Url::parse(&cli.alignment_bam).unwrap();
            IndexedReader::from_url(&url).unwrap()
        } else {
            IndexedReader::from_path(&cli.alignment_bam).unwrap()
        };
        
        info!("Successfully opened BAM file");

        bam
    }

pub fn process_bam_file(cli: &Cli, bam: &mut IndexedReader, chromosome: &str, windows: Vec<(usize, usize)>) -> Vec<Vec<fastq::Record>> {
    // Extract reads from the specified region
    let mut reads_list = Vec::new();
    for (start, end) in windows {
        let (reads, read_coordinates) = util::extract_reads(
            bam,
            &chromosome,
            start as u64,
            end as u64,
            Some(cli.quality), // min_mapping_quality
        ).unwrap();
        
        info!("Extracted {} reads from region", reads.len());

        // Extract read sequences from BAM file using utility function
        let read_sequence_dictionary = util::extract_read_sequences_from_region(
             bam,
            &chromosome,
            start as u64,
            end as u64,
        ).unwrap();

        for r in reads.iter() {
            let r_name = r.id().to_string().split('|').next().unwrap().to_string();
            let r_seq = read_sequence_dictionary.get(&r_name).unwrap();
            let r_start = read_coordinates.get(&r_name).unwrap().0;
            let r_end = read_coordinates.get(&r_name).unwrap().1 + 1;
            let reconstructed_seq = String::from_utf8_lossy(r.seq()).to_string();
            
            assert_eq!(r_seq[r_start as usize..r_end as usize], reconstructed_seq);
        }
        reads_list.push(reads);
    }
    reads_list
}

pub fn process_fasta_file(cli: &Cli, chromosome: &str, windows: Vec<(usize, usize)>) -> Vec<fastq::Record> {
    use bio::io::fasta::Reader as FastaReader;
    
    info!("Opening FASTA file: {}", cli.reference_fa);
    
    // Open FASTA file
    let mut fasta_reader = FastaReader::from_file(&cli.reference_fa)
        .expect("Failed to open FASTA file");
    
    let mut reference_seqs = Vec::new();
    
    // Read all sequences from FASTA
    for result in fasta_reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        let sequence = String::from_utf8_lossy(record.seq()).to_string();
        
        // info!("Found sequence: {} (length: {})", seq_id, sequence.len());
        
        // Check if this sequence matches the target chromosome
        if seq_id == chromosome {

            for (start, end) in windows.iter() {
                info!("Extracting region {}:{}-{} from chromosome {}", 
                    chromosome, start, end, seq_id);
                // Extract the specified region
                if start < &sequence.len() && end <= &sequence.len() {
                    let region_seq = sequence[*start..*end].to_string();
                    let record_id = format!("{}:{}-{}|reference|{}", 
                                        chromosome, start, end, cli.sampleid);
                    
                    let fastq_record = fastq::Record::with_attrs(
                        &record_id,
                        None,
                        region_seq.as_bytes(),
                        vec![30; region_seq.len()].as_slice() // Default quality score
                    );
                    
                    reference_seqs.push(fastq_record);
                    info!("Extracted reference sequence: {} bp", region_seq.len());
                } else {
                    warn!("Region coordinates out of bounds for sequence length {}", sequence.len());
                }
            }
        }
    }
    
    if reference_seqs.is_empty() {
        warn!("No matching chromosome '{}' found in FASTA file", chromosome);
    }
    
    info!("Extracted {} reference sequences", reference_seqs.len());
    reference_seqs
}


pub fn write_fasta_output(final_hap: HashMap<String, (String, Vec<String>, f64)>, output_path: &PathBuf) -> Result<()> {

    
    let mut file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;
    
    let mut index:i32 = 1;

    for (sequence, (cigar, read_names, allele_frequency)) in final_hap {
        // Use the first read name as the sequence ID, or create a hash-based ID
        let sequence_id = if !read_names.is_empty() {
            read_names.join(",")
        } else {
            "unknown".to_string()   
        };
        writeln!(file, ">Haplotype {}\tLength{}\tCigar {}\tSupportReads {}\tAlleleFrequency {:.2}\t{}", index, sequence.len(), cigar, read_names.len(), allele_frequency, sequence_id)?;
        writeln!(file, "{}", sequence)?;
        index += 1;
    }
    
    Ok(())
}

pub fn collapse_haplotypes(cli: &Cli, reads: &Vec<fastq::Record>, reference: &fastq::Record) -> Result<HashMap<String, (String, Vec<String>, f64)>> {
        // Read the reads records (name and sequence) into a vector.
        let mut read_dictionary: HashMap<String, Vec<String>> = HashMap::new();
        for read in reads{
            let contig_name = read.id().to_string();
            let contig = String::from_utf8_lossy(read.seq()).to_string();
            read_dictionary.entry(contig).or_default().push(contig_name);
        }
    
        let reference_seq = String::from_utf8_lossy(reference.seq()).to_string();
        let mut final_hap: HashMap<String, (String, Vec<String>, f64)> = HashMap::new();

        for (hap, vec) in read_dictionary.iter() {
            let allele_frequency = vec.len() as f64 / reads.len() as f64;
            if (vec.len() >= cli.min_reads as usize) && (allele_frequency >= cli.frequency_min as f64){
                let cigar = util::gap_open_aligner(&reference_seq, &hap);
                final_hap.insert(hap.clone(), (cigar, vec.clone(), allele_frequency));
            }
        }
        info!("Extracted {} haplotypes from region, total reads: {}", final_hap.len(), read_dictionary.len());
        Ok(final_hap)
    }

pub fn extract_haplotypes (cli: &Cli, chromosome: &str, windows: &Vec<(usize, usize)>) -> Result<(Vec<HashMap<String, (String, Vec<String>, f64)>>)> {
    // Validate inputs
    let mut bam = open_bam_file(cli);
    let reads_list:Vec<Vec<fastq::Record>> = process_bam_file(cli, &mut bam, chromosome, windows.clone());  
    let reference:Vec<fastq::Record> = process_fasta_file(cli, chromosome, windows.clone());
    let mut final_hap_list: Vec<HashMap<String, (String, Vec<String>, f64)>> = Vec::new();
    for (index, window) in windows.iter().enumerate() {
        let reads = &reads_list[index];
        let reference = &reference[index];
        let final_hap = collapse_haplotypes(cli, reads, reference)?;
        final_hap_list.push(final_hap);
    }

    Ok(final_hap_list)
}

pub fn get_node_edge_info(cli: &Cli, windows: &Vec<(usize, usize)>, final_hap_list: &Vec<HashMap<String, (String,Vec<String>, f64)>>) -> (HashMap<String, String>, HashMap<(String,String), String>) {
    let (chromosome, start, end) = util::split_locus(cli.locus.clone());
    // let gfa_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplograph.gfa", cli.output, cli.sampleid, chromosome, start, end));
    let mut node_info = HashMap::new();
    let mut edge_info = HashMap::new();
    for (index, window) in windows.iter().enumerate() {
        for (i, (final_haplotype_seq,(cigar, read_vector, allele_frequency))) in final_hap_list[index].iter().enumerate(){
            let haplotype_id = format!("Haplotype_{}|{}:{}-{}", i, chromosome, window.0, window.1);
            // let cigar = cigar_dict_list[index].get(final_haplotype_seq).unwrap();
            let read_vector_len = read_vector.len();
            node_info.insert(haplotype_id.clone(), format!("L\t{}\t{}\t{}\t{:.2}\t{}", final_haplotype_seq, cigar, read_vector_len, allele_frequency, read_vector.join(",")));

            // add edge information
            if index < windows.len() - 1 {
                let next_window = windows[index + 1];
                for (j, (next_final_haplotype_seq, (next_cigar,next_read_vector, next_allele_frequency))) in final_hap_list[index + 1].iter().enumerate(){
                    let next_haplotype_id = format!("Haplotype_{}|{}:{}-{}", j, chromosome, next_window.0, next_window.1);
                    // let next_cigar = cigar_dict_list[index + 1].get(next_final_haplotype_seq).unwrap();
                    let next_read_vector_len = next_read_vector.len();
                    let overlapping_reads = util::find_overlapping_reads(read_vector, next_read_vector);
                    let overlap_ratio = overlapping_reads.len() as f64 / (read_vector_len as f64).min( next_read_vector_len as f64);
                    // println!("readset1: {}, readset2: {}, overlap_ratio: {}", read_vector.len(), next_read_vector.len(), overlap_ratio);
                    if overlap_ratio > cli.frequency_min && overlapping_reads.len() > cli.min_reads as usize {
                        edge_info.insert((haplotype_id.clone(), next_haplotype_id.clone()), format!("E\t{}\t{}\tOverlapRatio:{:.3}\tOverlappingReads:{}", 
                            read_vector_len, next_read_vector_len, overlap_ratio, overlapping_reads.join(",")));
                    }
                    
                }
            }
        }
    }
    (node_info, edge_info)
}

pub fn write_gfa_output(
    node_file: &HashMap<String, String>,
    edge_info: &HashMap<(String, String), String>,
    output_filename: &PathBuf,
) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;
    writeln!(file, "H\tVN:Z:1.0")?;

    let mut node_output = Vec::new();
    for (haplotype_id, node_info) in node_file.iter() {
        let pos = haplotype_id.split("|").next().unwrap();
        let parts: Vec<&str> = node_info.split("\t").collect();
        if parts.len() >= 5 {
            let haplotype_seq = parts[1];
            let haplotype_cigar = parts[2];
            let read_num = parts[3];
            let allele_frequency = parts[4];
            
            let mut node_info_clone = BTreeMap::new();
            node_info_clone.insert("pos".to_string(), pos);
            // node_info_clone.insert("seq".to_string(), haplotype_seq);
            node_info_clone.insert("cigar".to_string(), haplotype_cigar);
            node_info_clone.insert("support_reads".to_string(), read_num);
            node_info_clone.insert("allele_frequency".to_string(), allele_frequency);
            // node_info_clone.insert("read_names".to_string(), read_names);
            // anchor_info_clone.seq = String::new();
            let json_string =
                serde_json::to_string(&node_info_clone).unwrap_or_else(|_| "{}".to_string());
            let formatted_string = format!("S\t{}\t{}\tPG:J:{}", haplotype_id, haplotype_seq, json_string);
            node_output.push(formatted_string);
        }

    }

    let mut link_output = Vec::new();

    for ((src, dst), edge_info_string) in edge_info.iter() {

        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", src, dst));


    }

    for s in node_output {
        writeln!(file, "{}", s)?;
    }
    for s in link_output {
        writeln!(file, "{}", s)?;
    }

    Ok(())
}

pub fn start(cli: &Cli, chromosome: &str, windows: &Vec<(usize, usize)>) -> Result<()> {
    let final_hap_list= extract_haplotypes(cli, chromosome, &windows)?;
    for (i, window) in windows.iter().enumerate() {
        let (start, end) = window;
        let haplotype_info = &final_hap_list[i];
        
        let final_hap_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplograph.fasta", cli.output, cli.sampleid, chromosome, start, end));
        write_fasta_output(haplotype_info.clone(), &final_hap_output)?;
    }

    info!("Haplotype reconstruction completed");
    Ok(())
}