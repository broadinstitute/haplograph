version 1.0

workflow Hifiasm {
    input {
        File whole_genome_bam
        File whole_genome_bai


        File reference_fa
        String prefix
        String locus

        Int hifiasm_mem
        Int hifiasm_thread
        Float min_freq
    }


    call SubsetBam {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            locus = locus,
            prefix = prefix
    }

    call hifiasm_asm {
        input:
            bam = SubsetBam.subsetbam,
            bai = SubsetBam.subsetbai,
            prefix = prefix + "_" + locus,
            num_cpus = hifiasm_thread,
            mem_gb = hifiasm_mem
    }


    output {

        File hifiasm_hap1 = hifiasm_asm.assembly_hap1
        File hifiasm_hap2 = hifiasm_asm.assembly_hap2

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

task SubsetBam {

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
        String prefix

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

