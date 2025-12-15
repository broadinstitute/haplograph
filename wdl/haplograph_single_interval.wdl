version 1.0

workflow Haplograph_single_interval {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
        File truth_vcf
        File truth_vcf_index
        String prefix
        String locus
        String truth_sample
        String sample_name
        Int coverage


    }



    call CalculateCoverage {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            locus = locus,
            prefix = prefix + "_" + locus
    }

    call haplograph {
        input:
            bam = CalculateCoverage.subsetbam,
            bai = CalculateCoverage.subsetbai,
            reference_fa = reference_fa,
            prefix = prefix,
            locus = locus
    }
    

    call Vcfdist {
        input:
            truth_sample = truth_sample,
            sample = prefix,
            genename = locus,
            coverage = coverage,
            eval_vcf = haplograph.germline_vcf_file,
            truth_vcf = truth_vcf,
            locus = locus,
            reference_fasta = reference_fa,
    }

    
    output {
        File precision_recall_summary = Vcfdist.precision_recall_summary_tsv
        File phasing_summary = Vcfdist.phasing_summary_tsv
    }
}

task haplograph {
    input {
        File bam
        File bai
        File reference_fa
        String prefix
        String locus
        Int windowsize
        Int minimal_supported_reads
        Int fold_threshold
        Float min_freq
        String extra_arg = ""
    }

    command <<<
        set -euxo pipefail
        /haplograph/target/release/haplograph haplograph -a ~{bam} \
                                                        -r ~{reference_fa} \
                                                        -s ~{prefix} \
                                                        -o ~{prefix} \
                                                        -l ~{locus} \
                                                        -v ~{min_freq} \
                                                        -m ~{minimal_supported_reads} \
                                                        -w ~{windowsize} \
                                                        -f gfa \
                                                        -c ~{fold_threshold}
                                                        ~{extra_arg}
        
        ls -l .
    >>>

    output {
        File graph_file = "~{prefix}.gfa"
        File asm_file = "~{prefix}.fasta"
        File germline_vcf_file = "~{prefix}.germline.vcf.gz"
        File somatic_vcf_file = "~{prefix}.somatic.vcf.gz"
        Array[File] methyl_bed = glob("*.bed")
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v3"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct DataTypeParameters {
    Int num_shards
    String map_preset
}

task CalculateCoverage {

    meta {
        description : "Subset a BAM file to a specified locus."
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
        samtools depth -r ~{locus} ~{prefix}.bam | awk '{sum+=$3} END {print sum/NR}' > coverage.txt

    >>>

    output {
        Float coverage = read_float("coverage.txt")
        File subsetbam =  "~{prefix}.bam"
        File subsetbai = " ~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task downsampleBam {

    input {
        File input_bam
        File input_bam_bai
        String basename
        Int desiredCoverage
        Float currentCoverage
        Int? preemptible_tries
    }

    meta {
        description: "Uses Picard to downsample to desired coverage based on provided estimate of coverage."
    }

    parameter_meta {
    }

    Float scalingFactor = desiredCoverage / currentCoverage


    command <<<
        set -eo pipefail
        if awk "BEGIN{ if((~{scalingFactor}) < 1.0) exit 0; else exit 1 }"; then
            gatk DownsampleSam -I ~{input_bam} -O ~{basename}_~{desiredCoverage}x.bam -R 7 -P ~{scalingFactor} -S ConstantMemory --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true
        else
            mv ~{input_bam} ~{basename}_~{desiredCoverage}x.bam
            mv ~{input_bam_bai} ~{basename}_~{desiredCoverage}x.bai
            echo "Total Coverage is lower than desiredCoverage"
        fi

    >>>
    runtime {
        preemptible: select_first([preemptible_tries, 5])
        memory: "8 GB"
        cpu: "2"
        disks: "local-disk 500 HDD"
        docker: "us.gcr.io/broad-gatk/gatk"
    }
    output {
        File downsampled_bam = "~{basename}_~{desiredCoverage}x.bam"
        File downsampled_bai = "~{basename}_~{desiredCoverage}x.bai"
    }
}

task split_into_bed {
    input {
        String locus
        Int bin_size
        Int pad_size
        String output_prefix

        Int? preemptible_tries
    }


    command <<<
        set -eo pipefail

        python - --locus ~{locus} \
                 --bin_size ~{bin_size} \
                 --pad_size ~{pad_size} \
                 --output_file ~{output_prefix} \
                 <<-'EOF'
        import gzip
        import argparse

        def split_locus(locus):
            chromosome, span = locus.split(":")
            start, end = span.split("-")
            return(chromosome, int(start), int(end))

        def split_locus_to_intervals(locus, bin_size=15000, pad_size=500):
            chromo, start, end = split_locus(locus)
            bin_num = (end - start)//bin_size
            intervals = [(chromo, start, start + bin_size + pad_size)]
            for i in range(1, bin_num):
                start_pos = start + i*bin_size - pad_size
                end_pos = start_pos + bin_size + pad_size
                intervals.append((chromo, start_pos, end_pos))
            intervals.append((chromo, start + bin_num*bin_size - pad_size, end))
            return(intervals)

        def write_bed_file(content, output_file):
            with open(output_file, "w") as f:
                for item in content:
                    l = "%s:%d-%d" % (item[0], item[1], item[2])
                    f.write(l+ "\n")

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--locus',
                                type=str)

            parser.add_argument('--output_file',
                                type=str)

            parser.add_argument('--bin_size',
                    type=int)

            parser.add_argument('--pad_size',
                    type=int)

            args = parser.parse_args()

            intervals = split_locus_to_intervals(args.locus, args.bin_size, args.pad_size)
            write_bed_file(intervals, args.output_file + ".txt")

        if __name__ == "__main__":
            main()
        EOF

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/slee/kage-lite:pr_29"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }

    output {
        Array[String] locuslist = read_lines("~{output_prefix}.txt")
    }
}

task ligate_vcfs {
    input {
        Array[File] vcfs
        String prefix
        Int disk_size_gb = 100
        Int mem_gb = 16
        Int cpu = 4
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # Set TMPDIR to current working directory to avoid /tmp write issues
        export TMPDIR=$PWD
        mkdir -p $TMPDIR

        for ff in ~{sep=' ' vcfs}; do bcftools index -t $ff; done

        bcftools concat \
            ~{sep=' ' vcfs} \
            --allow-overlaps \
            --ligate-force \
            --remove-duplicates \
            -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>
    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_tbi = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vcfdist:v1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

}

task Vcfdist {
    input {
        String truth_sample
        String sample
        String genename
        String coverage
        File eval_vcf
        File truth_vcf
        String locus
        File reference_fasta
        String? extra_args
        Int verbosity = 1

        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 32
        Int cpu = 8
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        bcftools index -t ~{truth_vcf}
        bcftools view -s ~{truth_sample} -r ~{locus} ~{truth_vcf} -Oz -o ~{sample}.~{locus}.base.vcf.gz
        bcftools index -t ~{sample}.~{locus}.base.vcf.gz

        bcftools index -t ~{eval_vcf}
        bcftools filter -e 'INFO/SOMATIC=1' ~{eval_vcf} -Oz -o ~{sample}.filtered.query.vcf.gz
        bcftools +fixploidy ~{sample}.filtered.query.vcf.gz -- -f 2 > ~{sample}.filtered.query.fixploidy.vcf
        bcftools view ~{sample}.filtered.query.fixploidy.vcf -Oz -o ~{sample}.filtered.query.fixploidy.vcf.gz
        bcftools index -t ~{sample}.filtered.query.fixploidy.vcf.gz

        vcfdist \
            ~{sample}.filtered.query.fixploidy.vcf.gz \
            ~{sample}.~{locus}.base.vcf.gz \
            ~{reference_fasta} \
            -v ~{verbosity} \
            ~{extra_args}

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.~{genename}.~{coverage}.$tsv; done
        mv summary.vcf ~{sample}.~{genename}.~{coverage}.summary.vcf
    >>>

    output {
        File summary_vcf = "~{sample}.~{genename}.~{coverage}.summary.vcf"
        File precision_recall_summary_tsv = "~{sample}.~{genename}.~{coverage}.precision-recall-summary.tsv"
        File query_tsv = "~{sample}.~{genename}.~{coverage}.query.tsv"
        File truth_tsv = "~{sample}.~{genename}.~{coverage}.truth.tsv"
        File phasing_summary_tsv = "~{sample}.~{genename}.~{coverage}.phasing-summary.tsv"
        File switchflips_tsv = "~{sample}.~{genename}.~{coverage}.switchflips.tsv"
        File phase_blocks_tsv = "~{sample}.~{genename}.~{coverage}.phase-blocks.tsv"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vcfdist:v1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}
