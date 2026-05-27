version 1.0

workflow Haplograph_pangenome {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File pan_fasta
        String prefix
        Int windowsize
    }

    call haplograph {
        input:
            bam = downsampleBam.downsampled_bam,
            bai = downsampleBam.downsampled_bai,
            reference_fa = reference_fa,
            prefix = prefix + "_" + gene_name + "_" + desiredCoverage,
            locus = locus,
            windowsize = windowsize
    }

    call haplograph_eval {
        input:
            truth_fasta = get_truth_haplotypes_from_annotation.fasta_file,
            query_fasta = haplograph.asm_file,
            prefix = prefix + "_" + gene_name + "_" + desiredCoverage,
    }


    output {
        Array[Float] bam_coverage = CalculateCoverage.coverage
        Array[Array[File]] gfa = haplograph.graph_file
        Array[Array[File]] fasta = haplograph.asm_file
        Array[Array[File]] vcf = haplograph.vcf_file
        Array[Array[File]] haplograph_eval_result = haplograph_eval.qv_scores
        Array[Array[VcfdistOutputs]] vcfdist_summary = VCFdist_germline.outputs
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

task pangenome_kmer_database {

    meta {
        description : "Construct kmer database from pangenome graph"
    }

    parameter_meta {
    }

    input {
        File pangenome_fa
        String prefix
        Int disk_size

        RuntimeAttr? runtime_attr_override
    }




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


