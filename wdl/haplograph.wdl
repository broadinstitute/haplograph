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
        Array[File] fasta = haplograph_asm.asm_file
    }
}

task haplograph_asm {
    input {
        File bam
        File bai
        File reference_fa
        String prefix
        String locus
    }

    command <<<
        set -euxo pipefail
        /haplograph/target/release/haplograph -a ~{bam} -r ~{reference_fa} -s ~{prefix} -l ~{locus} 
    >>>

    output {
        Array[File] asm_file = glob("*.fasta")
        
    }

    runtime {
        docker: "hangsuunc/haplograph:v1"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 300 SSD"
    }
}
