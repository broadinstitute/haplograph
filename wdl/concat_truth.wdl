version 1.0

workflow ConcatVCF {
    input {
        Array[File] vcfs
        Array[File] tbis
        String prefix
        String output_dir


    }
    call bcftools_concat {
        input:
            vcfs = vcfs,
            tbis = tbis,
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
        Array[File] tbis
        String prefix
        String output_dir
    }

    command <<<
        set -euxo pipefail

        bcftools concat \
            ~{sep=" " vcfs} \
            --allow-overlaps --remove-duplicates \
            -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz

        gsutil cp ~{prefix}.vcf.gz ~{output_dir}
        gsutil cp ~{prefix}.vcf.gz.tbi ~{output_dir}
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
