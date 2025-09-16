version 1.0

workflow Download {
    input {
        String address
        String output_address
        String prefix
    }

    call Download {
        input:
            address = address,
            output_bucket = output_address,
            prefix = prefix
    }


    output {
    }
}


task Download {
    input {
        String address
        String output_bucket
        String prefix
    }

    command <<<
        set -euxo pipefail
        wget ~{address}
        FILE_NAME=$(basename ~{address})
        gsutil cp $FILE_NAME ~{output_bucket}

        
    >>>

    output {
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 300 SSD"
    }
}
