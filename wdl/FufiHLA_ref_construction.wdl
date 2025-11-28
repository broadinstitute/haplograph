version 1.0
workflow FuFiHLA_ref_construction {
    input {
        String gcs_output_dir

    }
    

    call FuFiHLA_construction {
        input:
            gcs_output_dir = gcs_output_dir
    }

    output {

    }
}

task FuFiHLA_construction {
    input {
        String gcs_output_dir
    }

    command <<<
        set -euxo pipefail

        fufihla-ref-prep
        ls ref_data
        gsutil cp ref_data/ref.gene.fa.gz ~{gcs_output_dir}

    >>>

    output {
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/lr-fufihla:latest"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 200 SSD"
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