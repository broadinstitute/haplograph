version 1.0

workflow Immuannot {
    input {
        File hap1_bam
        File hap1_bai
        File hap2_bam
        File hap2_bai
        Array[String] locus_list
        String prefix
        String? billing_project
    }

    call ExtractAsm as hap1 {
        input:
            bam = hap1_bam,
            bai = hap1_bai,
            locus_list = locus_list,
            billing_project = billing_project,
            prefix = prefix + ".hap1"
    }

    call ExtractAsm as hap2 {
        input:
            bam = hap2_bam,
            bai = hap2_bai,
            locus_list = locus_list,
            billing_project = billing_project,
            prefix = prefix + ".hap2"
    }

    call immunoAnnotate as annot1 {
        input:
            target_asm = hap1.fasta,
            prefix = prefix + ".hap1"
    }

    call immunoAnnotate as annot2 {
        input:
            target_asm = hap2.fasta,
            prefix = prefix + ".hap2"
    }

    output {
        Array[File] hap1_annotation = annot1.annotation
        Array[File] hap2_annotation = annot2.annotation
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

task ExtractAsm {

    meta {
        description: "Subset a BAM file to specified loci and emit FASTA."
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai: "index for bam file"
        locus_list: "genomic loci to select"
        prefix: "prefix for output fasta file name"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        Array[String] locus_list
        String? billing_project
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4 * ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        export GCS_REQUESTER_PAYS_PROJECT="~{billing_project}"

        i=0
        for locus in ~{sep=" " locus_list}; do
            # Locus strings often contain ':' / '-'; use a safe per-interval filename.
            samtools view -bhX ~{bam} ~{bai} "$locus" | samtools fasta - > "locus_${i}.fa"
            i=$((i + 1))
        done

        cat locus_*.fa > ~{prefix}.fasta
    >>>

    output {
        File fasta = "~{prefix}.fasta"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
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

task immunoAnnotate {
    input {
        File target_asm
        String prefix
        Int threads = 2
    }

    command <<<
        set -euxo pipefail

        mkdir -p contigs output

        # Split multi-contig FASTA into one file per contig
        awk '/^>/ {
            if (f) close(f)
            name = substr($1, 2)
            gsub(/[^A-Za-z0-9_.-]/, "_", name)
            f = "contigs/" name ".fa"
            print > f
            next
        }
        { print > f }' ~{target_asm}

        shopt -s nullglob
        for fa in contigs/*.fa; do
            name=$(basename "${fa}" .fa)
            bash /Immuannot/scripts.pub.v3/immuannot.sh \
                -c "${fa}" \
                -r /Immuannot/Data-2024Feb02 \
                -o "output/~{prefix}.${name}" \
                -t ~{threads}
        done
    >>>

    output {
        Array[File] annotation = glob("output/*.gtf.gz")
    }

    Int disk_size = 100 + ceil(2 * size(target_asm, "GiB"))

    runtime {
        cpu: 4
        memory: "16 GiB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/immuannot:v1"
    }
}
