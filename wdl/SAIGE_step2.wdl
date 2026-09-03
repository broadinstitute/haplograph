version 1.0

workflow SAIGE_step2 {
    meta {
        description: "SAIGE step 2 only: single-variant and optional gene-set tests from pre-fit null models."
    }

    input {
        Array[File] null_model_rdas
        Array[File] variance_ratio_txts
        Array[String] phenotype_list
        File homoplasmic_vcf
        File homoplasmic_vcf_csi
        File? heteroplasmic_vcf
        File? heteroplasmic_vcf_csi
        File? GroupFile
        # Only supply these if Step 1 was run with --useSparseGRMtoFitNULL=TRUE. Step 1 writes a
        # "sparse"-tagged row into the varianceRatio file only in that case; if it is absent and
        # these are set, step2_SPAtests.R aborts with "sparse GRM is specified but the variance
        # ratio for sparse GRM was not estimated in Step 1".
        File? sparseGRM
        File? sparseGRM_IDlist
        String chromosome
        String trait_type
        String output_prefix
        String saige_docker
        String vcffield
        String memory = "8G"
        String? extra_argument
        Float single_variant_minimal_af = 0.01
        Int single_variant_min_mac = 20
        Float gene_set_minimal_af = 0
        Float gene_set_min_mac = 0.5

        # MAC category bounds for categorical variance ratios. These MUST match the values Step 1
        # was run with (--isCateVarianceRatio=TRUE), otherwise SAIGE silently applies each ratio to
        # the wrong MAC bin. SAIGE requires length(min_exclude) == length(max_include) + 1, and
        # length(min_exclude) == the number of "null" rows in the varianceRatio file.
        # Ignored by SAIGE when Step 1 used --isCateVarianceRatio=FALSE (single variance ratio).
        String cate_var_ratio_min_mac_exclude = "10,20.5"
        String cate_var_ratio_max_mac_include = "20.5"

        Int cpu = 2
        String disk_size = "local-disk 50 HDD"
        Int preemptible = 1
    }

    # Step 2a: single-variant test (homoplasmic VCF) for all phenotypes
    call RunStep2_singlevariant {
        input:
            vcf              = homoplasmic_vcf,
            vcf_csi          = homoplasmic_vcf_csi,
            variance_ratios  = variance_ratio_txts,
            GMMATmodelFiles  = null_model_rdas,
            sparseGRM        = sparseGRM,
            sparseGRM_IDlist = sparseGRM_IDlist,
            vcffield         = vcffield,
            chromo           = chromosome,
            phecode_list     = phenotype_list,
            trait_type       = trait_type,
            output_prefix    = output_prefix,
            minimal_af       = single_variant_minimal_af,
            min_mac          = single_variant_min_mac,
            cate_var_ratio_min_mac_exclude = cate_var_ratio_min_mac_exclude,
            cate_var_ratio_max_mac_include = cate_var_ratio_max_mac_include,
            memory           = memory,
            saige_docker     = saige_docker,
            cpu              = cpu,
            disk             = disk_size,
            preemptible      = preemptible,
            extra_argument   = extra_argument
    }

    # Step 2b: gene-set test (heteroplasmic VCF) for all phenotypes
    if (defined(heteroplasmic_vcf)){
        call RunStep2_geneset {
            input:
                vcf              = select_first([heteroplasmic_vcf, ""]),
                vcf_csi          = select_first([heteroplasmic_vcf_csi, ""]),
                variance_ratios  = variance_ratio_txts,
                GMMATmodelFiles  = null_model_rdas,
                vcffield         = vcffield,
                chromo           = chromosome,
                phecode_list     = phenotype_list,
                output_prefix    = output_prefix,
                minimal_af       = gene_set_minimal_af,
                min_mac          = gene_set_min_mac,
                GroupFile        = select_first([GroupFile, ""]),
                cate_var_ratio_min_mac_exclude = cate_var_ratio_min_mac_exclude,
                cate_var_ratio_max_mac_include = cate_var_ratio_max_mac_include,
                memory           = memory,
                saige_docker     = saige_docker,
                cpu              = cpu,
                disk             = disk_size,
                preemptible      = preemptible
        }
    }


    output {
        Array[File] singlevariant_outputs  = RunStep2_singlevariant.outputfiles
        Array[File]? geneset_outputs        = RunStep2_geneset.outputfiles
    }
}


# ---------------------------------------------------------------------------
# Task: fit null GLMM for every phecode in phecode_list
# ---------------------------------------------------------------------------
task RunFitNullGLMM {

    input {
        File plink_bed_file
        File plink_bim_file
        File plink_fam_file
        File sparseGRM
        File sparseGRM_IDlist
        File phenotype_file
        Array[String] phecode_list
        String trait_type
        String output_prefix
        String saige_docker

        # Must match the values passed to step2_SPAtests.R. See the workflow-level inputs.
        String cate_var_ratio_min_mac_exclude = "10,20.5"
        String cate_var_ratio_max_mac_include = "20.5"

        # Runtime parameters
        String memory
        Int cpu
        String disk
        Int preemptible
    }

    String inv_norm_arg = if (trait_type == "quantitative") then "--invNormalize=TRUE" else ""

    command <<<
        set -uxo pipefail

        phecodes=(~{sep=' ' phecode_list})
        successful_phecodes=()

        for phecode in "${phecodes[@]}"; do
            if step1_fitNULLGLMM.R \
                --bedFile="~{plink_bed_file}" \
                --bimFile="~{plink_bim_file}" \
                --famFile="~{plink_fam_file}" \
                --useSparseGRMtoFitNULL=TRUE \
                --sparseGRMFile="~{sparseGRM}" \
                --sparseGRMSampleIDFile="~{sparseGRM_IDlist}" \
                --phenoFile="~{phenotype_file}" \
                --phenoCol="${phecode}" \
                --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,sex,age,age2,age_sex,age2_sex \
                --qCovarColList=sex \
                --sampleIDColinphenoFile=person_id \
                --traitType=~{trait_type} \
                --isCateVarianceRatio=TRUE \
                --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
                --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
                --IsOverwriteVarianceRatioFile=TRUE \
                --outputPrefix="~{output_prefix}_step1Out_${phecode}" \
                ~{inv_norm_arg}; then
                successful_phecodes+=("${phecode}")
            else
                echo "WARNING: phecode ${phecode} failed, skipping" >&2
            fi
        done

        printf '%s\n' "${successful_phecodes[@]}" > successful_phecodes.txt
    >>>

    output {
        File        successful_phecode_list = "successful_phecodes.txt"
        Array[File] null_model_rdas         = glob("~{output_prefix}_step1Out_*.rda")
        Array[File] variance_ratio_txts     = glob("~{output_prefix}_step1Out_*.varianceRatio.txt")
    }

    runtime {
        docker:      saige_docker
        memory:      memory
        cpu:         cpu
        disks:       disk
        preemptible: preemptible
    }
}


# ---------------------------------------------------------------------------
# Task: step 2 single-variant SPA test for every phecode in phecode_list
# ---------------------------------------------------------------------------
task RunStep2_singlevariant {

    input {
        File vcf
        File vcf_csi
        Array[File] variance_ratios
        Array[File] GMMATmodelFiles
        # Set only when Step 1 used --useSparseGRMtoFitNULL=TRUE (see workflow-level inputs).
        File? sparseGRM
        File? sparseGRM_IDlist
        String vcffield
        String chromo
        Array[String] phecode_list
        String trait_type = "binary"
        String output_prefix
        Float minimal_af
        Int min_mac
        String cate_var_ratio_min_mac_exclude = "10,20.5"
        String cate_var_ratio_max_mac_include = "20.5"
        String? extra_argument

        # Runtime parameters
        String memory
        String saige_docker
        Int cpu
        String disk
        Int preemptible
    }

    command <<<
        set -euxo pipefail

        phecodes=(~{sep=' ' phecode_list})
        variance_ratios=(~{sep=' ' variance_ratios})
        model_files=(~{sep=' ' GMMATmodelFiles})

        for i in "${!phecodes[@]}"; do
            phecode="${phecodes[$i]}"
            variance_ratio="${variance_ratios[$i]}"
            model_file="${model_files[$i]}"

            step2_SPAtests.R \
                --vcfFile=~{vcf} \
                --vcfFileIndex=~{vcf_csi} \
                --vcfField=~{vcffield} \
                --AlleleOrder=ref-first \
                --chrom=~{chromo} \
                --SAIGEOutputFile=~{output_prefix}_step2Out_${phecode}_singlevar \
                --minMAF=~{minimal_af} \
                --minMAC=~{min_mac} \
                ~{"--sparseGRMFile=" + sparseGRM} \
                ~{"--sparseGRMSampleIDFile=" + sparseGRM_IDlist} \
                --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
                --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
                --GMMATmodelFile="${model_file}" \
                --varianceRatioFile="${variance_ratio}" \
                --is_Firth_beta=TRUE \
                --is_overrideflip=TRUE \
                --LOCO=FALSE \
                --is_fastTest=FALSE \
                --is_output_moreDetails=TRUE \
                ~{extra_argument}
        done
    >>>

    output {
        Array[File] outputfiles = glob("~{output_prefix}_step2Out_*_singlevar*")
    }

    runtime {
        docker:      saige_docker
        memory:      memory
        cpu:         cpu
        disks:       disk
        preemptible: preemptible
    }
}


# ---------------------------------------------------------------------------
# Task: step 2 gene-set test for every phecode in phecode_list
# ---------------------------------------------------------------------------
task RunStep2_geneset {
    meta {
        description: "Heteroplasmic variant association (all heteroplasmies) using SAIGE gene-set test — batch mode."
    }

    input {
        File vcf
        File vcf_csi
        Array[File] variance_ratios
        Array[File] GMMATmodelFiles
        String vcffield
        String chromo
        Array[String] phecode_list
        String output_prefix
        Float minimal_af
        Float min_mac
        File GroupFile
        String cate_var_ratio_min_mac_exclude = "10,20.5"
        String cate_var_ratio_max_mac_include = "20.5"

        # Runtime parameters
        String memory
        String saige_docker
        Int cpu
        String disk
        Int preemptible
    }

    command <<<
        set -euxo pipefail

        phecodes=(~{sep=' ' phecode_list})
        variance_ratios=(~{sep=' ' variance_ratios})
        model_files=(~{sep=' ' GMMATmodelFiles})

        for i in "${!phecodes[@]}"; do
            phecode="${phecodes[$i]}"
            variance_ratio="${variance_ratios[$i]}"
            model_file="${model_files[$i]}"

            step2_SPAtests.R \
                --vcfFile=~{vcf} \
                --vcfFileIndex=~{vcf_csi} \
                --vcfField=~{vcffield} \
                --AlleleOrder=ref-first \
                --chrom=~{chromo} \
                --SAIGEOutputFile=~{output_prefix}_step2Out_${phecode}_geneset \
                --minMAF=~{minimal_af} \
                --minMAC=~{min_mac} \
                --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
                --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
                --GMMATmodelFile="${model_file}" \
                --varianceRatioFile="${variance_ratio}" \
                --is_Firth_beta=TRUE \
                --LOCO=FALSE \
                --groupFile="~{GroupFile}" \
                --annotation_in_groupTest="non_coding_transcript_exon,start_lost,stop_gained;start_lost,missense,synonymous,frameshift,stop_gained,stop_lost,stop_retained,dloop" \
                --maxMAF_in_groupTest=0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5 \
                --pCutoffforFirth=0.05 \
                --is_output_markerList_in_groupTest=TRUE \
                --is_output_moreDetails=TRUE
        done
    >>>

    output {
        Array[File] outputfiles = glob("~{output_prefix}_step2Out_*_geneset*")
    }

    runtime {
        docker:      saige_docker
        memory:      memory
        cpu:         cpu
        disks:       disk
        preemptible: preemptible
    }
}
