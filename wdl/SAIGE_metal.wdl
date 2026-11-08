version 1.0

workflow SAIGE_metal {
    meta {
        description: "Fixed-effect inverse-variance meta-analysis of SAIGE_phewas_Step2 single-variant runs across 2 to 5 ancestries (METAL SCHEME STDERR). Pass each ancestry's singlevariant_outputs as study1_singlevar through study5_singlevar. study1 and study2 are required; study3, study4, and study5 are optional and must be filled in order, with no gaps. Phecodes are matched by the step2 file name. A phecode is pooled when at least two ancestries have results; chromosomes within an ancestry are concatenated. SAIGE BETA (the effect of Allele2) is pooled with weight 1/SE^2. META_BETA in the result is that estimate, on Allele2. Direction follows the ancestries included for that phecode, in study1 through study5 order; the studies column names those labels. Default labels are AFR, AMR, EUR, EAS, and SAS. The ancestries must be disjoint; SCHEME STDERR does not correct for sample overlap."
    }

    input {
        # RunStep2_singlevariant.outputfiles from each ancestry.
        # Supply 2 to 5 populations, in Direction order. study1 and study2 are
        # required; add study3, then study4, then study5. Do not skip a slot.
        Array[File] study1_singlevar
        Array[File] study2_singlevar
        Array[File]? study3_singlevar
        Array[File]? study4_singlevar
        Array[File]? study5_singlevar
        String study1_label = "AFR"
        String study2_label = "AMR"
        String study3_label = "EUR"
        String study4_label = "EAS"
        String study5_label = "SAS"
        # haplograph/scripts/saige_metal.py
        File saige_metal_py
        # Image from haplograph/docker/metal/Dockerfile. Needs python3 and metal on PATH.
        String metal_docker

        # Keep a marker only when every ancestry included for that phecode
        # contributed. Default keeps partial markers; n_studies counts how many did.
        Boolean both_studies_only = false

        String index_memory = "4G"
        String index_disk = "local-disk 200 HDD"
        String metal_memory = "8G"
        String metal_disk = "local-disk 100 HDD"
        Int cpu = 1
        Int preemptible = 1
    }

    call IndexPhecodes {
        input:
            study1_singlevar = study1_singlevar,
            study2_singlevar = study2_singlevar,
            study3_singlevar = study3_singlevar,
            study4_singlevar = study4_singlevar,
            study5_singlevar = study5_singlevar,
            study1_label = study1_label,
            study2_label = study2_label,
            study3_label = study3_label,
            study4_label = study4_label,
            study5_label = study5_label,
            saige_metal_py = saige_metal_py,
            metal_docker = metal_docker,
            memory = index_memory,
            cpu = cpu,
            disk = index_disk,
            preemptible = preemptible
    }

    scatter (bundle in IndexPhecodes.bundles) {
        call MetaPhecode {
            input:
                bundle = bundle,
                both_studies_only = both_studies_only,
                saige_metal_py = saige_metal_py,
                metal_docker = metal_docker,
                memory = metal_memory,
                cpu = cpu,
                disk = metal_disk,
                preemptible = preemptible
        }
    }

    output {
        File pairing_report = IndexPhecodes.pairing_report
        File direction_order = IndexPhecodes.direction_order
        File skipped_files = IndexPhecodes.skipped_files
        Array[File] meta_results = MetaPhecode.result
        Array[File] meta_info = MetaPhecode.info
        Array[File] meta_logs = MetaPhecode.log
        Array[File] meta_scripts = MetaPhecode.script
        Array[File] meta_qc = MetaPhecode.qc
    }
}


task IndexPhecodes {
    meta {
        description: "Match step2 single-variant files by phecode across the supplied ancestries (2 to 5)."
    }

    input {
        Array[File] study1_singlevar
        Array[File] study2_singlevar
        Array[File]? study3_singlevar
        Array[File]? study4_singlevar
        Array[File]? study5_singlevar
        String study1_label
        String study2_label
        String study3_label
        String study4_label
        String study5_label
        File saige_metal_py
        String metal_docker
        String memory
        Int cpu
        String disk
        Int preemptible
    }

    command <<<
        set -euo pipefail
        if [ "~{defined(study4_singlevar)}" = "true" ] && [ "~{defined(study3_singlevar)}" != "true" ]; then
            echo "study3_singlevar is required when study4_singlevar is set" >&2
            exit 1
        fi
        if [ "~{defined(study5_singlevar)}" = "true" ] && [ "~{defined(study4_singlevar)}" != "true" ]; then
            echo "study4_singlevar is required when study5_singlevar is set" >&2
            exit 1
        fi
        mkdir -p bundles
        cat > study1.files << 'EOF'
~{sep='\n' study1_singlevar}
EOF
        cat > study2.files << 'EOF'
~{sep='\n' study2_singlevar}
EOF
        cat > study3.files << 'EOF'
~{sep='\n' select_first([study3_singlevar, []])}
EOF
        cat > study4.files << 'EOF'
~{sep='\n' select_first([study4_singlevar, []])}
EOF
        cat > study5.files << 'EOF'
~{sep='\n' select_first([study5_singlevar, []])}
EOF
        # Full flag names. --study is an ambiguous prefix of --study1-list.
        args=(
            index
            --study1-list study1.files
            --study1-label "~{study1_label}"
            --study2-list study2.files
            --study2-label "~{study2_label}"
        )
        if [ "~{defined(study3_singlevar)}" = "true" ]; then
            args+=(--study3-list study3.files --study3-label "~{study3_label}")
        fi
        if [ "~{defined(study4_singlevar)}" = "true" ]; then
            args+=(--study4-list study4.files --study4-label "~{study4_label}")
        fi
        if [ "~{defined(study5_singlevar)}" = "true" ]; then
            args+=(--study5-list study5.files --study5-label "~{study5_label}")
        fi
        args+=(
            --bundle-dir bundles
            --pairing pairing.tsv
            --direction-order direction_order.tsv
            --skipped skipped_files.txt
        )
        python3 "~{saige_metal_py}" "${args[@]}"
    >>>

    output {
        Array[File] bundles = glob("bundles/*.tar.gz")
        File pairing_report = "pairing.tsv"
        File direction_order = "direction_order.tsv"
        File skipped_files = "skipped_files.txt"
    }

    runtime {
        docker: metal_docker
        memory: memory
        cpu: cpu
        disks: disk
        preemptible: preemptible
    }
}


task MetaPhecode {
    meta {
        description: "Concatenate chromosomes for one phecode and run METAL SCHEME STDERR."
    }

    input {
        File bundle
        Boolean both_studies_only
        File saige_metal_py
        String metal_docker
        String memory
        Int cpu
        String disk
        Int preemptible
    }

    command <<<
        set -euo pipefail
        python3 "~{saige_metal_py}" meta \
            --bundle "~{bundle}" \
            ~{true="--both-studies-only" false="" both_studies_only} \
            --out-dir .
    >>>

    output {
        File result = "metal_result.tsv"
        File info = "metal_result.info"
        File log = "metal.log"
        File script = "metal_script.txt"
        File qc = "metal_qc.tsv"
    }

    runtime {
        docker: metal_docker
        memory: memory
        cpu: cpu
        disks: disk
        preemptible: preemptible
    }
}
