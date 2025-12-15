version 1.0

workflow ConcatVCF {
    input {
        File VcfFile_list
        File TbiFile_list
        String input_dir
        String prefix
        String output_dir


    }
    Array[String] vcfs = read_lines(VcfFile_list)
    Array[String] tbis = read_lines(TbiFile_list)

    call bcftools_concat {
        input:
            vcfs = vcfs,
            vcf_tbis = tbis,
            prefix = prefix,
            output_dir = output_dir
    }

    output {
        File vcf = bcftools_concat.vcf
        File tbi = bcftools_concat.tbi
    }
}

task bcftools_concat {
    input {
        Array[File] vcfs
        String prefix
        String output_dir
    }

    command <<<
        set -euxo pipefail

        # Set TMPDIR to current working directory to avoid /tmp write issues
        export TMPDIR=$PWD
        mkdir -p $TMPDIR

        # Index all input VCF files
        for vcf in ~{sep=' ' vcfs}; do
            bcftools index "$vcf"
        done

        bcftools concat \
            ~{sep=" " vcfs} \
            --allow-overlaps --remove-duplicates \
            -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz

        bcftools sort ~{prefix}.vcf.gz -Oz -o ~{prefix}.sorted.vcf.gz
        bcftools index -t ~{prefix}.sorted.vcf.gz

        gsutil cp ~{prefix}.sorted.vcf.gz ~{output_dir}
        gsutil cp ~{prefix}.sorted.vcf.gz.tbi ~{output_dir}
    >>>

    output {
        File vcf = "~{prefix}.sorted.vcf.gz"
        File tbi = "~{prefix}.sorted.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 1000 SSD"
    }
}
