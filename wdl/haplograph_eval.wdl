version 1.0

workflow Haplograph_eval {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File truth_asm_1
        File truth_annotation_1
        File truth_asm_2
        File truth_annotation_2
        String gene_name
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

    call get_truth_haplotypes_from_annotation {
        input:
            truth_asm_1 = truth_asm_1,
            truth_asm_annotation_1 = truth_annotation_1,
            truth_asm_2 = truth_asm_2,
            truth_asm_annotation_2 = truth_annotation_2,
            reference_fa = reference_fa,
            gene_name = gene_name,
            prefix = prefix
    }

    call haplograph_eval {
        input:
            truth_fasta = get_truth_haplotypes_from_annotation.fasta_file,
            query_fasta = haplograph_asm.asm_file,
            seq_number = 2,
            prefix = prefix + locus
    }


    output {
        File gfa = haplograph_asm.graph_file
        File fasta = haplograph_asm.asm_file
        File eval = haplograph_eval.qv_scores
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
        File graph_file = "~{prefix}.gfa"
        File asm_file = "~{prefix}.fasta"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograp:v1"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}


task haplograph_eval {
    input {
        File truth_fasta
        File query_fasta
        Int seq_number
        String prefix
    }

    command <<<
        set -euxo pipefail

        /haplograph/target/release/haplograph evaluate -t ~{truth_fasta} \
                                                        -q ~{query_fasta} \
                                                        -s ~{seq_number} \
                                                        -o ~{prefix}.tsv
        
    >>>

    output {
        File qv_scores = "~{prefix}.tsv"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograp:v1"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}

task get_truth_haplotypes_from_annotation {
    input {
        File truth_asm_1
        File truth_asm_annotation_1
        File truth_asm_2
        File truth_asm_annotation_2
        File reference_fa
        String gene_name
        String prefix
    }

    command <<<
        set -euxo pipefail
        gunzip -c ~{truth_asm_1} > truth.hap1.fasta
        gunzip -c ~{truth_asm_2} > truth.hap2.fasta
        python - --truth_fasta_1 truth.hap1.fasta \
                 --truth_annotation_1 ~{truth_asm_annotation_1} \
                 --truth_fasta_2 truth.hap2.fasta \
                 --truth_annotation_2 ~{truth_asm_annotation_2} \
                 --gene_name ~{gene_name} \
                 --output_filename ~{prefix}.~{gene_name}.truth.fasta \
                 --output_header ~{gene_name} \
                 <<-'EOF'
        import pysam
        import gzip
        import argparse

        def get_hla_region(annotation_file, name):
            intervals = []
            with gzip.open(annotation_file, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    itemlist = line.strip().split("\t")
                    if itemlist[2] == "gene":
                        contig = itemlist[0]
                        start = int(itemlist[3])
                        end = int(itemlist[4])
                        for item in itemlist[8].split(";"):
                            if "gene_name" in item:
                                gene_name = item.split(" ")[2]
                                if name in gene_name:
                                    intervals.append((contig, start, end))
            return intervals  

        def get_truth_seq(fastafile, annotation_file, gene_name):
            fasta_records = pysam.Fastafile(fastafile)
            chromosomes = fasta_records.references
            truth_seq = []
            intervals = get_hla_region(annotation_file, gene_name)
            for chromosome in chromosomes:
                seq = fasta_records.fetch(chromosome)
                for (contig_name, start, end) in intervals:
                    if contig_name == chromosome:
                        truth_seq.append(seq[start:end])
            return truth_seq  


        def write_fasta(fasta_file, header, truth_seq):
            with open(fasta_file, "w") as f:
                for i, t_seq in enumerate(truth_seq):
                    f.write(">%s_Hap_%d\n" % (header, i+1))
                    char_length = 60
                    line_number = len(t_seq) // char_length

                    for i in range(line_number):
                        f.write(t_seq[i*char_length:(i+1)*char_length] + "\n")
                    f.write(t_seq[line_number*char_length:] + "\n")



        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--truth_fasta_1',
                                type=str)

            parser.add_argument('--truth_annotation_1',
                                type=str)

            parser.add_argument('--truth_fasta_2',
                                type=str)

            parser.add_argument('--truth_annotation_2',
                                type=str)

            parser.add_argument('--gene_name',
                                type=str)

            parser.add_argument('--output_filename',
                                type=str)
            parser.add_argument('--output_header',
                            type=str)

            args = parser.parse_args()

            truth_seq_1 = get_truth_seq(args.truth_fasta_1, args.truth_annotation_1, args.gene_name)
            truth_seq_2 = get_truth_seq(args.truth_fasta_2, args.truth_annotation_2, args.gene_name)
            truth_seq = truth_seq_1 + truth_seq_2
            write_fasta(args.output_filename, args.output_header, truth_seq)

        if __name__ == "__main__":
            main()
        EOF
        
    >>>

    output {
        File fasta_file = "~{prefix}.~{gene_name}.truth.fasta"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/slee/kage-lite:pr_29"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 100 SSD"
    }
}


task get_truth_haplotypes {
    input {
        File truth_hap1_bam
        File truth_hap1_bai
        File truth_hap2_bam
        File truth_hap2_bai
        File reference_fa
        String sampleid
        String locus
        String prefix
        String extra_arg
        Int windowsize
    }

    command <<<
        set -euxo pipefail

        /haplograph/target/release/haplograph haplotype -a ~{truth_hap1_bam} \
                                                -r ~{reference_fa} \
                                                -s ~{sampleid} \
                                                -o ~{prefix}.truth1 \
                                                -l ~{locus} \
                                                -w ~{windowsize} \
                                                -m 0 \
                                                -f 0 \
                                                -d fasta \
                                                ~{extra_arg}

        /haplograph/target/release/haplograph haplotype -a ~{truth_hap2_bam} \
                                                -r ~{reference_fa} \
                                                -s ~{sampleid} \
                                                -o ~{prefix}.truth2 \
                                                -l ~{locus} \
                                                -w ~{windowsize} \
                                                -m 0 \
                                                -f 0 \
                                                -d fasta \
                                                ~{extra_arg}
        cat *.fasta > ~{prefix}.truth.fasta
        
    >>>

    output {
        File fasta_file = "~{prefix}.truth.fasta"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograp:v1"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}
