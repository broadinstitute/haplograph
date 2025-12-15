version 1.0

workflow Haplograph_parameter_optimaize {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File truth_vcf
        File truth_vcf_tbi
        File reference_fa
        File reference_fai
        String prefix
        String locus
        String gene_name
        Int window
        Array[Int] coverages
        String truth_sample_name
        Float min_freq
    }

    call CalculateCoverage {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            locus = locus,
            prefix = prefix + "_" + gene_name
    }


    scatter (coverage in coverages) {
        
        call downsampleBam {input:
            input_bam = CalculateCoverage.subsetbam,
            input_bam_bai = CalculateCoverage.subsetbai,
            basename = prefix + "_" + gene_name,
            desiredCoverage = coverage,
            currentCoverage = CalculateCoverage.coverage,
            preemptible_tries = 0
        }

        call haplograph {
            input:
                bam = downsampleBam.downsampled_bam,
                bai = downsampleBam.downsampled_bai,
                reference_fa = reference_fa,
                prefix = prefix + "_" + gene_name + "_" + window,
                locus = locus,
                windowsize = window,
                min_freq = min_freq
        }

        call VCFEval as VCFEVAL {
            input:
                query_vcf = haplograph.somatic_vcf_file,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                base_vcf = truth_vcf,
                base_vcf_index = truth_vcf_tbi,
                query_output_sample_name = prefix + "_" + coverage + "_" + gene_name,
                prefix = prefix + "_" + coverage + "_" + gene_name,
                truth_sample = truth_sample_name,
                locus = locus
        }

    }
    


    output {
        Array[File] gfa = haplograph.graph_file
        Array[File] fasta = haplograph.asm_file
        Array[File] germline_vcf = haplograph.germline_vcf_file
        Array[File] somatic_vcf = haplograph.somatic_vcf_file
        Array[File] vcfeval_summary = VCFEVAL.summary_statistics
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


task haplograph_eval {
    input {
        File truth_fasta
        File query_fasta
        String prefix
    }

    command <<<
        set -euxo pipefail
        /haplograph/target/release/haplograph evaluate -t ~{truth_fasta} \
                                                        -q ~{query_fasta} \
                                                        -s 2 \
                                                        -o ~{prefix}.tsv
        
    >>>

    output {
        File qv_scores = "~{prefix}.tsv"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v3"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}


task get_truth_haplotypes {
    input {
        File truth_hap1_bam
        File truth_hap1_bai
        File truth_hap2_bam
        File truth_hap2_bai
        String locus
        String prefix
    }

    command <<<
        set -euxo pipefail

        /haplograph/target/release/haplograph extract -b ~{truth_hap1_bam} \
                                                -o ~{prefix}.truth1 \
                                                -l ~{locus} 

        /haplograph/target/release/haplograph extract -b ~{truth_hap2_bam} \
                                                -o ~{prefix}.truth2 \
                                                -l ~{locus}

        cat ~{prefix}.truth1.fasta ~{prefix}.truth2.fasta > ~{prefix}.truth.fasta
        
    >>>

    output {
        File fasta_file = "~{prefix}.truth.fasta"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v3"
        memory: "4 GB"
        cpu: 1
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
        cpu_cores:          1,
        mem_gb:             10,
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

task parseBed {

    input {
        File bed
        String output_prefix

        Int? preemptible_tries
    }


    command <<<
        set -eo pipefail

        python - --bed_file ~{bed} \
                 --output_file ~{output_prefix} \
                 <<-'EOF'
        import gzip
        import argparse


        def import_bed(bed_file):
            locus_list = []
            gene_list = []
            if bed_file.endswith("gz"):
                with gzip.open(bed_file, "r") as f:
                    for line in f:
                        itemlist = line.strip().split("\t")
                        locus = "%s:%d-%d" % (itemlist[0], int(itemlist[1]), int(itemlist[2]))
                        locus_list.append(locus)
                        gene_list.append(itemlist[3])
            else:
                with open(bed_file, "r") as f:
                    for line in f:
                        itemlist = line.strip().split("\t")
                        locus = "%s:%d-%d" % (itemlist[0], int(itemlist[1]), int(itemlist[2]))
                        locus_list.append(locus)
                        gene_list.append(itemlist[3])
            return (locus_list, gene_list)


        def write_file(content, output_file):
            with open(output_file, "w") as f:
                for item in content:
                    f.write(f"{item}\n")

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--bed_file',
                                type=str)

            parser.add_argument('--output_file',
                                type=str)

            args = parser.parse_args()

            locus_list, gene_list = import_bed(args.bed_file)
            write_file(locus_list, args.output_file + "_locus.txt")
            write_file(gene_list, args.output_file + "_gene.txt")


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
        Array[String] locuslist = read_lines("~{output_prefix}_locus.txt")
        Array[String] genelist= read_lines("~{output_prefix}_gene.txt")
    }
}

struct VcfdistOutputs {
    File summary_vcf
    File precision_recall_summary_tsv
    File precision_recall_tsv
    File query_tsv
    File truth_tsv
    File phasing_summary_tsv
    File switchflips_tsv
    File superclusters_tsv
    File phase_blocks_tsv
}

struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
}

task Vcfdist {
    input {
        String truth_sample
        String sample
        String genename
        String coverage
        File eval_vcf
        File? eval_vcf_index
        File truth_vcf
        File? truth_vcf_index
        String locus
        File reference_fasta
        String? extra_args
        Int verbosity = 1

        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail
        bcftools index -t ~{truth_vcf}
        bcftools view -s ~{truth_sample} -r ~{locus} ~{truth_vcf} -Oz -o ~{sample}.~{locus}.base.vcf.gz
        bcftools index -t ~{sample}.~{locus}.base.vcf.gz

        bcftools index -t ~{eval_vcf}

        vcfdist \
            ~{eval_vcf} \
            ~{sample}.~{locus}.base.vcf.gz \
            ~{reference_fasta} \
            -v ~{verbosity} \
            ~{extra_args}

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.~{genename}.~{coverage}.$tsv; done
        mv summary.vcf ~{sample}.~{genename}.~{coverage}.summary.vcf
    >>>

    output {
        VcfdistOutputs outputs = {
            "summary_vcf": "~{sample}.~{genename}.~{coverage}.summary.vcf",
            "precision_recall_summary_tsv": "~{sample}.~{genename}.~{coverage}.precision-recall-summary.tsv",
            "precision_recall_tsv": "~{sample}.~{genename}.~{coverage}.precision-recall.tsv",
            "query_tsv": "~{sample}.~{genename}.~{coverage}.query.tsv",
            "truth_tsv": "~{sample}.~{genename}.~{coverage}.truth.tsv",
            "phasing_summary_tsv": "~{sample}.~{genename}.~{coverage}.phasing-summary.tsv",
            "switchflips_tsv": "~{sample}.~{genename}.~{coverage}.switchflips.tsv",
            "superclusters_tsv": "~{sample}.~{genename}.~{coverage}.superclusters.tsv",
            "phase_blocks_tsv": "~{sample}.~{genename}.~{coverage}.phase-blocks.tsv"
        }
        File base_vcf = "~{sample}.~{locus}.base.vcf.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vcfdist:v1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

task hifiasm_asm{
    input{
        File bam
        File bai
        String prefix
        Int num_cpus
        Int mem_gb
    }

    Int disk_size = 10 + ceil(2 * size(bam, "GiB"))

    command <<<

        set -euxo pipefail

        samtools fastq ~{bam} > ~{prefix}.fastq

        /truvari/hifiasm-0.25.0/hifiasm -o ~{prefix} -t 4 ~{prefix}.fastq
        awk '/^S/{print ">"$2;print $3}' ~{prefix}.bp.hap1.p_ctg.gfa > ~{prefix}.bp.hap1.p_ctg.fa
        awk '/^S/{print ">"$2;print $3}' ~{prefix}.bp.hap2.p_ctg.gfa > ~{prefix}.bp.hap2.p_ctg.fa

        cat ~{prefix}.bp.hap1.p_ctg.fa ~{prefix}.bp.hap2.p_ctg.fa > ~{prefix}.hifiasm.fa
    >>>

    output{
        File assembly_hap1="~{prefix}.bp.hap1.p_ctg.fa"
        File assembly_hap2="~{prefix}.bp.hap2.p_ctg.fa"
        File asm_file = "~{prefix}.hifiasm.fa"
    }
    runtime {
        cpu: num_cpus
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/hifiasm:0.25.0"
    }    
}

task VCFEval {
    input {
        # Input VCF Files
        File query_vcf
        File reference_fa
        File reference_fai
        File base_vcf
        File base_vcf_index
        String query_output_sample_name
        String truth_sample
        String prefix
        String locus

        # Filtering params
        Int? min_indel_length = 50  # Exclude indels smaller than this length (default 50bp for HiFi)

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB") + 2 * size(base_vcf, "GB") + size(reference_fa, "GB")) + 50,
                                                  "cpu": 8, "memory": 16}
    }

    command <<<
        set -xeuo pipefail
        bcftools view -s ~{truth_sample} -r ~{locus} ~{base_vcf} -Oz -o ~{prefix}.~{locus}.base.vcf.gz
        bcftools index -t ~{prefix}.~{locus}.base.vcf.gz

        # Compress and Index vcf files
        bcftools view ~{query_vcf} -O z -o ~{query_output_sample_name}.vcf.gz
        bcftools index -t ~{query_output_sample_name}.vcf.gz

        # 
        
        # rtg vcfeval
        rtg format -o rtg_ref ~{reference_fa}
        rtg vcfeval \
            -b ~{prefix}.~{locus}.base.vcf.gz  \
            -c ~{query_output_sample_name}.vcf.gz \
            -o reg \
            -t rtg_ref \
            --squash-ploidy \
            --sample ALT,ALT 

        mkdir output_dir
        cp reg/summary.txt output_dir/~{query_output_sample_name}_summary.txt
        
    
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1-tmp"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + " GB"
    }

    output {
        File summary_statistics = "output_dir/~{query_output_sample_name}_summary.txt"
        
    }
}

task aardvark {
    input {
        File truth_vcf
        File truth_vcf_index
        File query_vcf
        File reference_fa
        String prefix
        String locus
        String truth_sample
        File bedfile
    }

    command <<<
        set -euxo pipefail

        bcftools index -t ~{query_vcf}
        
        bcftools view -s ~{truth_sample} -r ~{locus} ~{truth_vcf} -Oz -o ~{prefix}.~{locus}.base.vcf.gz
        bcftools index -t ~{prefix}.~{locus}.base.vcf.gz

        /truvari/aardvark-v0.9.0-x86_64-unknown-linux-gnu/aardvark compare --min-variant-gap 1000 \
                                                        --reference ~{reference_fa} \
                                                        --truth-vcf ~{prefix}.~{locus}.base.vcf.gz  \
                                                        --query-vcf ~{query_vcf} \
                                                        --regions ~{bedfile} \
                                                        --output-dir ~{prefix}

        cp ~{prefix}/summary.tsv ~{prefix}.~{locus}.summary.tsv
        cp ~{prefix}/truth.vcf.gz ~{prefix}.~{locus}.truth.vcf.gz
        cp ~{prefix}/query.vcf.gz ~{prefix}.~{locus}.query.vcf.gz
        
    >>>

    output {
        File summary_file = "~{prefix}.~{locus}.summary.tsv"
        File output_truth_vcf = "~{prefix}.~{locus}.truth.vcf.gz"
        File output_query_vcf = "~{prefix}.~{locus}.query.vcf.gz"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/aardvark:0.9.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task concat_vcfs {
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
        bcftools sort ~{prefix}.vcf.gz -Oz -o ~{prefix}.sorted.vcf.gz
        bcftools index -t ~{prefix}.sorted.vcf.gz
    >>>
    output {
        File merged_vcf = "~{prefix}.sorted.vcf.gz"
        File merged_vcf_tbi = "~{prefix}.sorted.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vcfdist:v1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

}