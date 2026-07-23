version 1.0

workflow Haplograph_eval_regular_genes {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File hap1_asm_bam
        File hap1_asm_bai
        File hap2_asm_bam
        File hap2_asm_bai
        File truth_vcf
        File? truth_vcf_tbi
        File reference_fa

        String sample
        String truth_sample = ""
        String output_prefix
        String coverage_tag = ""
        Int maximal_locus_size = 200000
        Int overlap_bp = 0
        Float min_mean_depth = 2.0
        Int window_size = 100
        Int vcfdist_verbosity = 1
        Int haplograph_threads = 4
        String extra_vcfdist_args = ""
        Boolean skip_haplograph = false
        Boolean verbose = false

        File gene_bed
        Array[Int] desiredCoverages
        Int hifiasm_mem
        Int hifiasm_thread
        Float min_freq
        String haplograph_docker = "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:dev"
    }



    call parseBed {
        input:
            bed = gene_bed,
            output_prefix = output_prefix
    }

    scatter (pair in zip(parseBed.locuslist, parseBed.genelist)) {
        String locus = pair.left
        String gene_name = pair.right

        call CalculateCoverage {
            input:
                bam = whole_genome_bam,
                bai = whole_genome_bai,
                locus = locus,
                prefix = output_prefix + "_" + gene_name
        }

        call parse_locus {
            input:
                locus = locus,
                maximal_locus_size = maximal_locus_size
        }

        call get_truth_haplotypes {
            input:
                truth_hap1_bam = hap1_asm_bam,
                truth_hap1_bai = hap1_asm_bai,
                truth_hap2_bam = hap2_asm_bam,
                truth_hap2_bai = hap2_asm_bai,
                locus = locus,
                prefix = output_prefix + "_" + gene_name
        }

        scatter (desiredCoverage in desiredCoverages) {
            
            call downsampleBam {input:
                input_bam = CalculateCoverage.subsetbam,
                input_bam_bai = CalculateCoverage.subsetbai,
                basename = output_prefix + "_" + gene_name,
                desiredCoverage = desiredCoverage,
                currentCoverage = CalculateCoverage.coverage,
                preemptible_tries = 0
            }

            call run_haplograph_benchmark as haplograph {
                input:
                    bam = downsampleBam.downsampled_bam,
                    bai = downsampleBam.downsampled_bai,
                    reference_fa = reference_fa,
                    sample = sample,
                    output_prefix = output_prefix+ "_" + gene_name + "_" + desiredCoverage,
                    locus = locus,
                    needs_merge = parse_locus.needs_merge,
                    window_size = window_size,
                    maximal_locus_size = maximal_locus_size,
                    overlap_bp = overlap_bp,
                    min_mean_depth = min_mean_depth,
                    verbose = verbose,
                    haplograph_docker = haplograph_docker,
                    thread = haplograph_threads
            }

            # Only run haplograph-dependent tasks when coverage was sufficient.
            if (!haplograph.low_coverage) {
                call haplograph_eval {
                    input:
                        truth_fasta = get_truth_haplotypes.fasta_file,
                        query_fasta = haplograph.asm_file,
                        prefix = output_prefix + "_" + gene_name + "_" + desiredCoverage,
                }

                call Vcfdist as VCFdist_germline {
                    input:
                        sample = sample,
                        eval_vcf = haplograph.germline_vcf,
                        truth_vcf = truth_vcf,
                        locus = locus,
                        reference_fasta = reference_fa,
                        coverage = desiredCoverage,
                        genename = gene_name,
                        extra_args = ""
                }
            }

            call hifiasm_asm {
                input:
                    bam = downsampleBam.downsampled_bam,
                    bai = downsampleBam.downsampled_bai,
                    prefix = output_prefix + "_" + gene_name + "_" + desiredCoverage,
                    num_cpus = hifiasm_thread,
                    mem_gb = hifiasm_mem
            }

            call haplograph_eval as hifiasm_eval {
                input:
                    truth_fasta = get_truth_haplotypes.fasta_file,
                    query_fasta = hifiasm_asm.asm_file,
                    prefix = output_prefix + "_" + gene_name + "_" + desiredCoverage,
            }
        }
    }


    output {
        Array[Float] bam_coverage = CalculateCoverage.coverage
        Array[Array[File?]] gfa = haplograph.graph_file
        Array[Array[File]] fasta = haplograph.asm_file
        Array[Array[File]] vcf = haplograph.germline_vcf
        Array[Array[File?]] haplograph_eval_result = haplograph_eval.qv_scores
        Array[Array[File]] hifiasm_eval_result = hifiasm_eval.qv_scores
        Array[Array[VcfdistOutputs?]] vcfdist_summary = VCFdist_germline.outputs
    }
}


task parse_locus {
    input {
        String locus
        Int maximal_locus_size
    }

    command <<<
        set -euo pipefail
        python3 - <<'PY'
        import sys

        locus = "~{locus}"
        maximal = int("~{maximal_locus_size}")

        _chrom, span = locus.split(":", 1)
        start_s, end_s = span.split("-", 1)
        start = int(start_s)
        end = int(end_s)
        locus_len = end - start
        if locus_len <= 0:
            sys.exit(f"invalid locus span {start}-{end} => {locus_len} bp")

        locus_tag = locus.replace(":", "_").replace("-", "_")
        needs_merge = locus_len > maximal

        with open("locus_len.txt", "w") as fh:
            fh.write(str(locus_len))
        with open("locus_tag.txt", "w") as fh:
            fh.write(locus_tag)
        with open("needs_merge.txt", "w") as fh:
            fh.write("true" if needs_merge else "false")
        PY
    >>>

    output {
        Int locus_len = read_int("locus_len.txt")
        String locus_tag = read_string("locus_tag.txt")
        Boolean needs_merge = read_boolean("needs_merge.txt")
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/slee/kage-lite:pr_29"
        memory: "1 GiB"
        cpu: 1
        disks: "local-disk 10 HDD"
    }
}



task run_haplograph_benchmark {
    meta {
        description: "Run haplograph haplograph; call haplograph dev-tools merge when the locus exceeds maximal_locus_size."
    }

    input {
        File bam
        File bai
        File reference_fa
        String sample
        String output_prefix
        String locus
        Boolean needs_merge
        Int window_size
        Int maximal_locus_size
        Int overlap_bp
        Int thread
        Float min_mean_depth
        Boolean verbose
        String haplograph_docker

        RuntimeAttr? runtime_attr_override
    }

    String verbose_flag = if (verbose) then "--verbose" else ""
    Int disk_gb = 20 + ceil(size([bam, reference_fa], "GiB"))

    command <<<
        set -uxo pipefail

        export TMPDIR="${PWD}"
        mkdir -p "${TMPDIR}"

        HAPLOGRAPH=/haplograph/target/release/haplograph

        LOW_COVERAGE=false
        ${HAPLOGRAPH} haplograph \
            --alignment-bam ~{bam} \
            --reference-fa ~{reference_fa} \
            --sampleid ~{sample} \
            --output-prefix ~{output_prefix} \
            --locus ~{locus} \
            --window-size ~{window_size} \
            --maximal-locus-size ~{maximal_locus_size} \
            --overlap-bp ~{overlap_bp} \
            --min-mean-depth ~{min_mean_depth} \
            --threads ~{thread} \
            ~{verbose_flag} \
            || { LOW_COVERAGE=true; }

        if [ "${LOW_COVERAGE}" = "true" ] || [ ! -f "~{output_prefix}.fasta" ]; then
            echo "Insufficient coverage or no output produced; skipping locus." >&2
            touch ~{output_prefix}.fasta
            touch ~{output_prefix}.germline.vcf.gz
            touch ~{output_prefix}.somatic.vcf.gz
            echo "true" > low_coverage.txt
        else
            echo "false" > low_coverage.txt
            if [ "~{needs_merge}" = "true" ]; then
                ${HAPLOGRAPH} dev-tools merge \
                    --output-prefix ~{output_prefix} \
                    --locus ~{locus} \
                    --reference-fa ~{reference_fa} \
                    --sampleid ~{sample} \
                    --maximal-locus-size ~{maximal_locus_size} \
                    --overlap-bp ~{overlap_bp} \
                    ~{verbose_flag}
            fi
        fi

    >>>

    output {
        File? graph_file = "~{output_prefix}.gfa"
        File asm_file = "~{output_prefix}.fasta"
        File germline_vcf = "~{output_prefix}.germline.vcf.gz"
        File? germline_vcf_tbi = "~{output_prefix}.germline.vcf.gz.tbi"
        File somatic_vcf_file = "~{output_prefix}.somatic.vcf.gz"
        File? somatic_vcf_tbi = "~{output_prefix}.somatic.vcf.gz.tbi"
        Array[File] methyl_bed = glob("~{output_prefix}.Hap.*.bed")
        Boolean low_coverage = read_boolean("low_coverage.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             haplograph_docker
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
        RuntimeAttr? runtime_attr_override
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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:dev"
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


task haplograph_eval {
    input {
        File truth_fasta
        File query_fasta
        String prefix
        RuntimeAttr? runtime_attr_override
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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:dev"
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


task get_truth_haplotypes {
    input {
        File truth_hap1_bam
        File truth_hap1_bai
        File truth_hap2_bam
        File truth_hap2_bai
        String locus
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        /haplograph/target/release/haplograph extract -b ~{truth_hap1_bam} \
                                                -l ~{locus} \
                                                -o ~{prefix}.truth1 \
                                                -s "~{prefix}_hap1"

        /haplograph/target/release/haplograph extract -b ~{truth_hap2_bam} \
                                                -l ~{locus} \
                                                -o ~{prefix}.truth2 \
                                                -s "~{prefix}_hap2"

        truth_files=()
        [[ -f ~{prefix}.truth1.fasta ]] && truth_files+=(~{prefix}.truth1.fasta)
        [[ -f ~{prefix}.truth2.fasta ]] && truth_files+=(~{prefix}.truth2.fasta)

        if [[ ${#truth_files[@]} -gt 0 ]]; then
            cat "${truth_files[@]}" > ~{prefix}.truth.fasta
        else
            touch ~{prefix}.truth.fasta
        fi

    >>>

    output {
        File fasta_file = "~{prefix}.truth.fasta"
        
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:dev"
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
        RuntimeAttr? runtime_attr_override
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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk"
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
        bcftools filter -e 'INFO/SOMATIC=1' ~{eval_vcf} -Oz -o ~{sample}.filtered.query.vcf.gz
        bcftools index -t ~{sample}.filtered.query.vcf.gz

        vcfdist \
            ~{sample}.filtered.query.vcf.gz \
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
        RuntimeAttr? runtime_attr_override
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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_gb,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/hifiasm:0.25.0"
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

