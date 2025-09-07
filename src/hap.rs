use log::{info, warn};
use anyhow::{Result as AnyhowResult, Context};
use bio::io::fastq;
use rust_htslib::faidx::Reader;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use std::path::{PathBuf};
use url::Url;
use std::collections::{HashMap, BTreeMap, HashSet};
use std::io::Write;
use std::fs::File;
use crate::Cli;
use crate::util;
use crate::intervals;
use std::error::Error;
use indicatif::{ProgressBar, ProgressStyle};

pub fn start(bam: &mut IndexedReader, chromosome: &str, windows: &Vec<(usize, usize)>, locus:String, reference_fa: &Vec<fastq::Record>, sampleid: &String, min_reads: usize, frequency_min: f64, primary_only: bool, output_prefix: &String) -> AnyhowResult<()> {

    let (chromosome, start, end) = util::split_locus(locus);
    for (i, window) in windows.iter().enumerate() {
        let (start, end) = window;
        let haplotype_info = intervals::start(bam, &reference_fa, &chromosome, *start, *end, &sampleid, min_reads as usize, frequency_min, primary_only, true).unwrap();
    }

    info!("Haplotype reconstruction completed");
    Ok(())
}