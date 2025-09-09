version 1.0

workflow Haplograph {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
        String prefix
        String locus


    }
    call haplograph_asm {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            reference_fa = reference_fa,
            prefix = prefix,
            locus = locus
    }


    output {
        File fasta = haplograph_asm.asm_file
    }
}

task haplograph_asm {
    input {
        File bam
        File bai
        File reference_fa
        String prefix
        String locus
        Int windowsize = 1000
        String extra_arg = ""
    }

    command <<<
        set -euxo pipefail
        /haplograph/target/release/haplograph haplotype -a ~{bam} \
                                                        -r ~{reference_fa} \
                                                        -s ~{prefix} \
                                                        -o ~{prefix} \
                                                        -l ~{locus} \
                                                        -w ~{windowsize} \
                                                        -d gfa \
                                                        ~{extra_arg}
        /haplograph/target/release/haplograph assemble -g ~{prefix}.gfa -o ~{prefix}.fasta
        
    >>>

    output {
        File asm_file = "~{prefix}.fasta"
        
    }

    runtime {
        docker: "hangsuunc/haplograph:v1"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}
