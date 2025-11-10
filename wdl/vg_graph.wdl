
version 1.0

workflow VG_extract_bam {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File gfa_file
        String prefix
        String locus
        Int giraffe_threads
        String sample_id
        String giraffe_extra_arg
        String graph_ref_path
    }

    call SubsetBam {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            locus = locus,
            prefix = prefix
    }

    call vg_graph_index {
        input:
            graph = gfa_file,
            prefix = prefix
    }

    call vg_graph_alignment {
        input:
            graph_index = vg_graph_index.graph_index,
            graph_dist = vg_graph_index.graph_dist,
            graph_min = vg_graph_index.graph_min,
            graph_gfa = gfa_file,
            graph_zipcode = vg_graph_index.graph_zipcode,
            input_fq = SubsetBam.subsetfq,
            thread_num = giraffe_threads,
            in_sample_name = sample_id,
            extra_args = giraffe_extra_arg,
            out_prefix = prefix + "_vg_align",
            ref_path = graph_ref_path

    }

    
    output {
        File graph_bam = vg_graph_alignment.bam_file
        File graph_bai = vg_graph_alignment.bam_index_file

    }
}

task SubsetBam {

    input {
        File bam
        File bai
        String locus
        String prefix
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam

        samtools fastq -F 0x900 -OT RG,BC,Mm,Ml -n ~{prefix}.bam > ~{prefix}.fastq
 
    >>>

    output {
        File subsetbam =  "~{prefix}.bam"
        File subsetbai = "~{prefix}.bam.bai"
        File subsetfq = "~{prefix}.fastq"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vg:1.69.0"
        cpu: 1
        memory: "4GB"
        disks: "local-disk " + 100 + " HDD"
        preemptible: 0
    }

}



task vg_graph_index{
    input{
        File graph
        String prefix
    }

    command <<<

    set -euxo pipefail
    /vg/vg autoindex --workflow lr-giraffe \
                     --prefix ~{prefix} \
                     --gfa ~{graph}

    >>>

    output {
        File graph_index= "~{prefix}.giraffe.gbz"
        File graph_dist = "~{prefix}.dist"
        File graph_min = "~{prefix}.longread.withzip.min"
        File graph_zipcode = "~{prefix}.longread.zipcodes"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vg:1.69.0"
        cpu: 1
        memory: "4GB"
        disks: "local-disk " + 100 + " HDD"
        preemptible: 0
    }
}

task vg_graph_alignment {

    input {
        File graph_index
        File graph_dist
        File graph_min
        File graph_gfa
        File graph_zipcode
        File input_fq
        Int thread_num
        Boolean input_is_gam = false
        String read_type = "hifi"
        String in_sample_name
        String extra_args
        String out_prefix
        String ref_path
    }
    Int half_cores = thread_num / 2

    command <<<
        set -euxo pipefail

        /vg/vg giraffe --parameter-preset ~{read_type} \
                       --progress --track-provenance \
                       -t ~{thread_num} \
                       -Z ~{graph_index} \
                       -d ~{graph_dist} \
                       -m ~{graph_min} \
                       -z ~{graph_zipcode} \
                       -f ~{input_fq} \
                       -o gam > output.gam

        /vg/vg surject \
                -n ~{ref_path} \
                -x ~{graph_index} \
                -t ~{thread_num} \
                --multimap \
                --bam-output ~{true="" false="--gaf-input" input_is_gam} \
                --sample ~{in_sample_name} \
                --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:hifi PU:unit1" \
                ~{extra_args} \
                output.gam | samtools sort --threads ~{half_cores} \
                                            -O BAM > ~{out_prefix}.bam
        samtools index -b ~{out_prefix}.bam
        
    >>>

    output {
        File bam_file = "~{out_prefix}.bam"
        File bam_index_file = "~{out_prefix}.bam.bai"
    }

    #########################

    runtime {
        docker:  "hangsuunc/vg:v1"
        cpu: 16
        memory: "64GB"
        disks: "local-disk " + 200 + " HDD"
        preemptible: 0
    }
    
}