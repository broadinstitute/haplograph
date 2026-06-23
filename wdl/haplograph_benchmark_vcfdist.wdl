version 1.0

workflow haplograph_benchmark_vcfdist {
    meta {
        description: "Localize BAM, run haplograph (+ merge for large loci), compare to truth with vcfdist."
    }

    input {
        File alignment_bam
        File alignment_bai
        File reference_fa
        String locus
        File truth_vcf
        File? truth_vcf_index
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
        File? query_vcf
        String haplograph_docker = "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:dev"
        String samtools_docker = "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
        String vcfdist_docker = "us.gcr.io/broad-dsp-lrma/hangsuunc/vcfdist:v1"
    }

    String effective_truth_sample = if (truth_sample == "") then sample else truth_sample

    call parse_locus {
        input:
            locus = locus,
            maximal_locus_size = maximal_locus_size
    }


    if (!skip_haplograph) {
        call CalculateCoverage {
            input:
                bam = alignment_bam,
                bai = alignment_bai,
                locus = locus,
                prefix = output_prefix + "_" + locus
        }

        call run_haplograph_benchmark {
            input:
                bam = CalculateCoverage.subsetbam,
                bai = CalculateCoverage.subsetbai,
                reference_fa = reference_fa,
                sample = sample,
                output_prefix = output_prefix,
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
    }

    File final_query_vcf = select_first([run_haplograph_benchmark.germline_vcf, query_vcf])

    call run_vcfdist {
        input:
            query_vcf = final_query_vcf,
            truth_vcf = truth_vcf,
            truth_vcf_index = truth_vcf_index,
            truth_sample = effective_truth_sample,
            sample = sample,
            locus = locus,
            locus_tag = parse_locus.locus_tag,
            coverage_tag = coverage_tag,
            reference_fa = reference_fa,
            verbosity = vcfdist_verbosity,
            extra_args = extra_vcfdist_args,
            docker = vcfdist_docker
    }

    output {
        Int locus_len = parse_locus.locus_len
        Boolean merged = parse_locus.needs_merge
        File query_germline_vcf = final_query_vcf
        File? query_germline_vcf_index = run_haplograph_benchmark.germline_vcf_tbi
        VcfdistOutputs vcfdist = run_vcfdist.outputs
        File precision_recall_summary = run_vcfdist.outputs.precision_recall_summary_tsv
        File phasing_summary = run_vcfdist.outputs.phasing_summary_tsv
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



task run_haplograph_benchmark {
    meta {
        description: "Run haplograph haplograph; call haplograph merge when the locus exceeds maximal_locus_size."
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
        set -euxo pipefail

        export TMPDIR="${PWD}"
        mkdir -p "${TMPDIR}"

        HAPLOGRAPH=/haplograph/target/release/haplograph

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
            ~{verbose_flag}

        if [ "~{needs_merge}" = "true" ]; then
            ${HAPLOGRAPH} merge \
                --output-prefix ~{output_prefix} \
                --locus ~{locus} \
                --reference-fa ~{reference_fa} \
                --sampleid ~{sample} \
                --maximal-locus-size ~{maximal_locus_size} \
                --overlap-bp ~{overlap_bp} \
                ~{verbose_flag}
        fi

    >>>

    output {
        File germline_vcf = "~{output_prefix}.germline.vcf.gz"
        File? germline_vcf_tbi = "~{output_prefix}.germline.vcf.gz.tbi"
        Array[File] segment_outputs = glob("~{output_prefix}_*")
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

task run_vcfdist {
    input {
        File query_vcf
        File truth_vcf
        File? truth_vcf_index
        String truth_sample
        String sample
        String locus
        String locus_tag
        String coverage_tag
        File reference_fa
        Int verbosity
        String extra_args
        String docker

        RuntimeAttr? runtime_attr_override
    }

    String out_tag = sample + "." + locus_tag + "." + coverage_tag
    Int disk_gb = 10 + ceil(size([query_vcf, truth_vcf, reference_fa], "GiB"))

    command <<<
        set -euxo pipefail

        export TMPDIR="${PWD}"
        mkdir -p "${TMPDIR}"

        QUERY_VCF="~{query_vcf}"
        bcftools index -t "${QUERY_VCF}"

        if [ ! -f "~{truth_vcf}.tbi" ] && [ ! -f "~{truth_vcf}.csi" ]; then
            bcftools index -t ~{truth_vcf}
        fi
        bcftools view -s ~{truth_sample} -r ~{locus} ~{truth_vcf} -Oz -o truth.subset.vcf.gz
        bcftools index -t truth.subset.vcf.gz

        bcftools filter -e 'INFO/SOMATIC=1' "${QUERY_VCF}" -Oz -o query.filtered.vcf.gz
        bcftools index -t query.filtered.vcf.gz

        vcfdist \
            query.filtered.vcf.gz \
            truth.subset.vcf.gz \
            ~{reference_fa} \
            -v ~{verbosity} \
            ~{extra_args}

        for tsv in *.tsv; do
            mv "${tsv}" "~{out_tag}.${tsv}"
        done
        mv summary.vcf "~{out_tag}.summary.vcf"
    >>>

    output {
        VcfdistOutputs outputs = {
            "summary_vcf": "~{out_tag}.summary.vcf",
            "precision_recall_summary_tsv": "~{out_tag}.precision-recall-summary.tsv",
            "precision_recall_tsv": "~{out_tag}.precision-recall.tsv",
            "query_tsv": "~{out_tag}.query.tsv",
            "truth_tsv": "~{out_tag}.truth.tsv",
            "phasing_summary_tsv": "~{out_tag}.phasing-summary.tsv",
            "switchflips_tsv": "~{out_tag}.switchflips.tsv",
            "superclusters_tsv": "~{out_tag}.superclusters.tsv",
            "phase_blocks_tsv": "~{out_tag}.phase-blocks.tsv"
        }
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             docker
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
