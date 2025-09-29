use log::{info, warn};
use anyhow::{Result as AnyhowResult, Context};
use bio::io::fastq;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use rust_htslib::bcf::{self, index as BcfIndex, record::GenotypeAllele};
use std::path::{PathBuf};
use std::collections::{HashMap, BTreeMap, HashSet};
use crate::Cli;
use crate::util;
use crate::intervals;
use indicatif::{ProgressBar, ProgressStyle};

pub fn get_all_haplotypes_to_vcf(bam: &mut IndexedReader, windows: &Vec<(String, usize, usize)>,  reference_fa: &Vec<fastq::Record>, sampleid: &String, min_reads: usize, frequency_min: f64, primary_only: bool, output_prefix: &String, default_file_format: &String) -> AnyhowResult<()> {
    for (i, window) in windows.iter().enumerate() {
        let (chromosome, start, end) = window;
        let haplotype_info = intervals::start(bam, &reference_fa, &chromosome, *start, *end, &sampleid, min_reads as usize, frequency_min, primary_only, false).unwrap();
    }
    Ok(())
}

pub fn start(bam: &mut IndexedReader, windows: &Vec<(String, usize, usize)>,  reference_fa: &Vec<fastq::Record>, sampleid: &String, min_reads: usize, frequency_min: f64, primary_only: bool, output_prefix: &String, default_file_format: &String) -> AnyhowResult<()> {
    if default_file_format == "fasta" { 
        for (i, window) in windows.iter().enumerate() {
            let (chromosome, start, end) = window;
            let haplotype_info = intervals::start(bam, &reference_fa, &chromosome, *start, *end, &sampleid, min_reads as usize, frequency_min, primary_only, true).unwrap();
        }
    } else if default_file_format == "vcf" {
        let mut header = bcf::Header::new();
        let mut reference_sequence = String::new();
        let mut reference_name = String::new();
        
        for record in reference_fa.iter() {
            let referencename = record.id().to_string();
            let referencelength = record.seq().len() as u64;
            header.push_record(format!("##contig=<ID={},length={}>\n", referencename, referencelength).as_bytes());
        }
        // header.push_record(format!("##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n").as_bytes());
        header.push_record(format!("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n").as_bytes());
        header.push_record(format!("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n").as_bytes());
        header.push_record(format!("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Alternative Allele Depth\">\n").as_bytes());
        header.push_record(format!("##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant Allele Frequency\">\n").as_bytes());
        header.push_sample(sampleid.as_bytes());

        let output_file = PathBuf::from(format!("{}.vcf.gz", output_prefix));
        let mut writer = bcf::Writer::from_path(output_file.clone(), &header, false, bcf::Format::Vcf).expect("Failed to create BCF writer");
    
        for (i, window) in windows.iter().enumerate() {
            let (chromosome,start, end) = window;
            let referene_record= reference_fa.iter().find(|r| r.id().to_string() == *chromosome).unwrap();
            let reference_sequence = String::from_utf8_lossy(referene_record.seq()).to_string();
            let reference_name = referene_record.id().to_string();
            let header_view = writer.header();
            let reference_id = header_view.name2rid(reference_name.as_bytes()).unwrap();

            let mut record = writer.empty_record();
            let reference_seq = reference_sequence[*start..*end].to_string();
            let haplotype_info = intervals::start(bam, &reference_fa, &chromosome, *start, *end, &sampleid, min_reads as usize, frequency_min, primary_only, false).unwrap();
            let mut record_list = Vec::new();
            record_list.push((&reference_seq, &0.0, 1 as usize));
            let mut coverage = 0;
            let mut haplotype_info_sorted = haplotype_info.iter().collect::<Vec<_>>();
            haplotype_info_sorted.sort_by(|a, b| a.0.cmp(&b.0));
            for (hap, (cigar, reads, allele_frequency)) in haplotype_info_sorted.iter() {
                coverage += reads.len();
                if **hap == reference_seq{
                    continue
                }
                record_list.push((hap, allele_frequency, reads.len()));
            }
            if record_list.len() <= 1 {
                continue
            }
            record.set_rid(Some(reference_id as u32));
            record.set_pos(*start as i64);
            record.set_id(b".");
            record.set_alleles(&record_list.iter().map(|v| v.0.as_bytes()).collect::<Vec<&[u8]>>()).expect("Failed to set alleles");
            record.push_format_integer(b"DP", &[coverage as i32]).expect("Failed to set DP format field");
            record.push_format_integer(b"AD", &record_list.iter().skip(1).map(|v| v.2 as i32).collect::<Vec<i32>>()).expect("Failed to set AD format field");
            record.push_format_float(b"VAF", &record_list.iter().skip(1).map(|v| *v.1 as f32).collect::<Vec<f32>>()).expect("Failed to set VAF format field");
            writer.write(&record).expect("Failed to write record");
        }
    }


    info!("Haplotype reconstruction completed");
    Ok(())
}