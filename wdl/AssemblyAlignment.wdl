version 1.0

workflow Minimap2_AlignAsm {
    input {
        File query_fasta
        File reference_fa
        String prefix
    }

    call AlignAsm {
        input:
            query_fa = query_fasta,
            reference_fa = reference_fa,
            prefix = prefix
    }


    output {
        File alignmentfile = AlignAsm.alignment_file
        File alignmentindex = AlignAsm.alignment_index
    }
}


task AlignAsm {
    input {
        File query_fa
        File reference_fa
        String prefix
    }

    command <<<
        set -euxo pipefail
        minimap2 -ayYL --MD --eqx -x asm5 -R '@RG\tID:as,\tSM:asm' -t 4 ~{reference_fa} ~{query_fa} -o ~{prefix}.sam
        samtools sort -O BAM --write-index -o ~{prefix}.bam ~{prefix}.sam
        
    >>>

    output {
        File alignment_file = "~{prefix}.bam"
        File alignment_index = "~{prefix}.bam.bai"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 300 SSD"
    }
}
