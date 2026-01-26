version 1.0

workflow MidasHLA_v1{
    meta{
        description: "a workflow for HLA association using MidasHLA"
    }
    input{
        File hla_allele_type
        File kir_allele_type
        File phecodex_table
        File meta_data
        String phecode
        Int hla_resolution


        String midashla_memory = "8G"
        Int midashla_cpu = 2
        String midashla_disk = "local-disk 50 HDD"
        Int midashla_preemptible = 1

    }

    call PreprocessFiles{input:
        hla_allele_file = hla_allele_type,
        kir_allele_file = kir_allele_type,
        phenotype_data = phecodex_table,
        meta_data = meta_data,
        phecode = phecode,
        disk_size = 100
    }

    call RunMidasHLA {input:
        hla_alleles_types = PreprocessFiles.hla_alleles,
        kir_alleles_types = PreprocessFiles.kir_alleles,
        phenotype_file = PreprocessFiles.phenotype,
        phecode = phecode,
        hla_resolution = hla_resolution,
        memory = midashla_memory,
        cpu = midashla_cpu,
        disk = midashla_disk,
        preemptible = midashla_preemptible
    }

    output {
        File phenotype = PreprocessFiles.phenotype
        Array[File] association_results = RunMidasHLA.outputfiles
        Array[File] figures = RunMidasHLA.outputfigures
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



task PreprocessFiles {
    input {
        File hla_allele_file
        File kir_allele_file
        File phenotype_data
        File meta_data
        String phecode
        Int disk_size

    }

    command <<<
        set -euxo pipefail

        python - --hla_allele ~{hla_allele_file} \
                --kir_allele ~{kir_allele_file} \
                --phenotype_data ~{phenotype_data} \
                --meta_data ~{meta_data} \
                --phecode ~{phecode} \
                <<-'EOF'
        import gzip
        import argparse
        import json
        import pandas as pd
        import numpy
        import random
        from datetime import date

        def get_phenotype_info(filename, pid_list, phecode, chunksize):
            df_phenotype = pd.DataFrame()
            with pd.read_csv(filename, chunksize=chunksize, dtype = str) as reader:
                for chunk in reader: 
                    df_chunk = chunk[chunk["person_id"].isin(pid_list)]
                    df_chunk_pheno = df_chunk[['person_id', phecode]]
                    df_phenotype = pd.concat([df_phenotype, df_chunk_pheno], axis = 0, ignore_index=True)
            df_phenotype = df_phenotype.drop_duplicates()
            return df_phenotype

        def get_demographic_data(metadata, pid_list, chunksize, columns = ['person_id',"sex_at_birth","ancestry_pred_other",'date_of_birth',"PC1", "PC2", "PC3", "PC4", "PC5"]):
            index = 0
            df_all = pd.DataFrame()
            with pd.read_csv(metadata, chunksize=chunksize, dtype = str) as reader:
                for chunk in reader:
                    df_chunk = chunk[chunk["person_id"].isin(pid_list)]
                    df_sub = df_chunk[columns]
                    df_all = pd.concat([df_all, df_sub], axis = 0, ignore_index = True)
            df_all = df_all.drop_duplicates() 
            today = date.today()
            df_age = pd.to_datetime(df_all['date_of_birth'])
            df_all['age'] = [int(today.year - t.year) for t in df_age]
            df_all['age_groups'] = [int((today.year - t.year)/10) * 10 for t in df_age]
            return df_all

        def get_case_pheno_info(df_pheno, df_demographic, phecode):
            case_sample = df_pheno[df_pheno[phecode] == "TRUE"]
            case_sample = case_sample.set_index(['person_id'])
            case_pids = case_sample.index.tolist()
            case_demographic = df_demographic[df_demographic['person_id'].isin(case_pids)]
            case_demographic = case_demographic.set_index(['person_id'])
            df_case_pheno = pd.concat([case_sample, case_demographic], axis = 1)
            return df_case_pheno

        def match_control_samples(df_pheno, df_demographic, df_case_pheno, phecode, columns = ["ancestry_pred_other", "sex_at_birth", 'age_groups']):
            control_sample = df_pheno[df_pheno[phecode] == "FALSE"]
            control_sample = control_sample.set_index(['person_id'])
            control_pid = control_sample.index.tolist()
            control_demographic = df_demographic[df_demographic['person_id'].isin(control_pid)]
            control_demographic = control_demographic.set_index(['person_id'])
            df_control_pheno = pd.concat([control_sample, control_demographic], axis = 1)

            df_case_info = df_case_pheno[columns].value_counts()
            match_control_sample = []
            total = 0
            for key, count in df_case_info.items():
                total += count
                anc,sex,age = key
                df_match_control_anc = df_control_pheno[df_control_pheno["ancestry_pred_other"] == anc]
                df_match_control_anc_sex = df_match_control_anc[df_match_control_anc['sex_at_birth'] == sex]
                df_match_control_anc_sex_age = df_match_control_anc_sex[df_match_control_anc_sex['age_groups'] == age]
                if len(df_match_control_anc_sex_age) >= count:
                    match_id_list = df_match_control_anc_sex_age.index.tolist()
                elif len(df_match_control_anc_sex) >= count:
                    match_id_list = df_match_control_anc_sex.index.tolist()
                    print(key, count, "anc, sex")
                elif len(df_match_control_anc) >= count:
                    match_id_list = df_match_control_anc.index.tolist()
                    print(key, count, "anc")
                else:
                    match_id_list = control_sample_info.index.tolist()
                    print(key, count, "no-match")
                match_control_sample += random.sample(match_id_list, k=count)
            df_control_pheno_matched = df_control_pheno.loc[match_control_sample, :]
            return df_control_pheno_matched

        def filter_inconsistent_calls(df_allele):
            for c in df_allele.columns:
                genename = c.split("_")[0]
                df_col = df_allele[c]
                for index in df_col.index:
                    v = df_col[index]
                    if pd.isna(v):
                        continue
                    if not v.startswith(genename):
                        df_allele.loc[index, c] = numpy.nan
                        print(genename, v, c, index)
            return df_allele

        def reformat_hla_alleles(filename, df_merged, phecode, geneset = {"A", "B", "C", "DR", "DB", "DQ", "DP"}):
            df_allele = pd.read_csv(filename, sep = ",", index_col=0)
            df_allele.index = df_allele.index.astype(str)
            experiment_sample = df_merged.index.tolist()
            columns = []
            for c_ in df_allele.columns:
                c = c_.split("_")[0]
                if len(c)>1:
                    c = c[:2]
                if c in geneset:
                    columns.append(c_)
            df_allele_experiment = df_allele.loc[experiment_sample, columns]
            for c in columns:
                df_allele_experiment[c] = df_allele_experiment[c].str.replace('HLA-', '', regex=False)
            df_allele_experiment_filtered = filter_inconsistent_calls(df_allele_experiment)
            df_allele_experiment_filtered = df_allele_experiment_filtered.dropna(axis=1, how='all')
            df_allele_experiment_filtered.to_csv("%s.hla.csv" % phecode, index=True)


        def reformat_kir_alleles(filename, df_merged, phecode):
            df_kir = pd.read_csv(filename, sep = ",", index_col=0)
            df_kir.index = df_kir.index.astype(str)
            experiment_sample = df_merged.index.astype(str).tolist()
            kir_dict = {}
            for index, row in df_kir.iterrows():
                for allele in row.values:
                    if pd.isna(allele):
                        continue
                    a = allele.split("*")[0]
                    kir_dict[index] = kir_dict.get(index, {})
                    kir_dict[index][a] = kir_dict[index].get(a, 0) + 1
            df_kir_allele = pd.DataFrame(kir_dict).T
            df_kir_allele_experiment = df_kir_allele.loc[experiment_sample,:]
            df_kir_allele_experiment = df_kir_allele_experiment.fillna(0)
            df_kir_allele_experiment.to_csv("%s.kir.csv" % phecode, index=True)


        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--hla_allele_file',
                                type=str)
            parser.add_argument('--kir_allele_file',
                                type=str)
            parser.add_argument('--phenotype_data',
                                type=str)
            parser.add_argument('--meta_data',
                                type=str)
            parser.add_argument('--phecode',
                                type=str)

            args = parser.parse_args()
            phecode = args.phecode

            # load HLA alleles
            df_hla = pd.read_csv(args.hla_allele_file, sep = ",", header=0, index_col=0)
            pid_list = df_hla.index.astype(str).tolist()

            # load KIR alleles
            df_kir = pd.read_csv(args.kir_allele_file, sep = ",", index_col=0)
            df_kir.index = df_kir.index.astype(str)

            # load Phenotypes
            df_pheno = get_phenotype_info(args.phenotype_data, pid_list, phecode, chunksize = 1000)
            df_demographic = get_demographic_data(args.meta_data, pid_list, chunksize = 1000)

            df_case_pheno = get_case_pheno_info(df_pheno, df_demographic, phecode)

            df_control_samples = match_control_samples(df_pheno, df_demographic, df_case_pheno, phecode, columns = ["ancestry_pred_other", "sex_at_birth", 'age_groups'])
            df_merged = pd.concat([df_case_pheno,df_control_samples])
            replacement_map = {'TRUE': 1, 'FALSE': 0}
            df_merged = df_merged.replace(replacement_map)
            df_merged.to_csv("%s.pheno.csv" % phecode, index=True)

            reformat_hla_alleles(args.hla_allele_file, df_merged, phecode)
            reformat_kir_alleles(args.kir_allele_file, df_merged, phecode)

        if __name__ == "__main__":
            main()

        EOF

    >>>

    output {
        File phenotype = "~{phecode}.pheno.csv"
        File hla_alleles = "~{phecode}.hla.csv"
        File kir_alleles = "~{phecode}.kir.csv"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/midashla:v0"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " SSD"
    }
}

task RunMidasHLA {

  input {
    File hla_alleles_types
    File kir_alleles_types
    File phenotype_file
    String phecode
    Int hla_resolution = 4

    # Runtime parameters
    String memory
    Int cpu
    String disk
    Int preemptible
  }
  command <<<
    set -euxo pipefail

    echo '
            args <- commandArgs(trailingOnly = TRUE)
            hla_file <- args[1]
            kir_file <- args[2]
            phenotype_file <- args[3]
            phecode <- args[4]
            hla_resolution <- as.integer(args[5])
            
            library(midasHLA)
            library(readr)
            library(tidyr)
            library(ggplot2)
            library(dplyr)
            library(stringr)
            library(ggrepel)

            # read hla alleles
            hla_data <- read_csv(hla_file)
            write.table(hla_data, paste(hla_file, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
            dat_HLA <- readHlaCalls( file = paste(hla_file, ".tsv"), resolution = hla_resolution, na.strings = c("None", "-", "NA") ) 
            
            # read phenotype
            pheno <- read.table(phenotype_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
            names(pheno)[1] <- "ID"

            # read kir alleles
            kir_data <- read_csv(kir_file)
            names(kir_data)[1] <-"ID"
            df_reordered <- kir_data[, c("ID", "KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A", "KIR2DS3", "KIR2DS5", "KIR2DS4", "KIR2DS1", "KIR3DL2")]
            names(df_reordered) = c("ID", "KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5", "KIR2DS3", "KIR2DS5", "KIR2DS4", "KIR2DS1", "KIR3DL2")
            write.table(df_reordered, paste(kir_file, ".tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
            kir_calls <- readKirCalls(paste(kir_file, ".tsv"))
            
            # determine covariates (sex specific phenotypes)
            if ("sex_at_birth" %in% names(pheno)) {
              sex_values <- pheno$sex_at_birth[!is.na(pheno$sex_at_birth)]
              unique_sex <- unique(sex_values)
              if (length(unique_sex) < 2) {
                cat("Warning: sex_at_birth has no variation (all values are:", paste(unique_sex, collapse = ", "), ")\n")
                cat("Removing sex_at_birth from the model formula to avoid GLM errors.\n")
                # Create a flag to exclude sex from formulas
                exclude_sex <- TRUE
              } else {
                cat("sex_at_birth has variation:", paste(unique_sex, collapse = ", "), "\n")
                exclude_sex <- FALSE
              }
            } else {
              cat("Warning: sex_at_birth column not found in phenotype data.\n")
              exclude_sex <- TRUE
            }

            if ("age" %in% names(pheno)) {
              age_values <- pheno$age[!is.na(pheno$age)]
              unique_age <- unique(age_values)
              if (length(unique_age) < 2) {
                cat("Warning: age has no variation (all values are:", paste(unique_age, collapse = ", "), ")\n")
                cat("Removing age from the model formula to avoid GLM errors.\n")
                # Create a flag to exclude sex from formulas
                exclude_age <- TRUE
              } else {
                cat("age has variation:", paste(unique_age, collapse = ", "), "\n")
                exclude_age <- FALSE
              }
            } else {
              cat("Warning: age column not found in phenotype data.\n")
              exclude_age <- TRUE
            }


            # Build covariate list (exclude sex_at_birth if no variation)
            covariates <- c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "term")
            if (!exclude_sex && !exclude_age) {
              covariates <- c("age", "sex_at_birth", "PC1", "PC2", "PC3", "PC4", "PC5", "term")
            }else if (exclude_sex) {
              covariates <- c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "term")
            }else if (exclude_age){
              covariates <- c("sex_at_birth","PC1", "PC2", "PC3", "PC4", "PC5", "term")
            }else {
              covariates <- c("PC1", "PC2", "PC3", "PC4", "PC5", "term")
            
            }
            print(covariates)


            # hla association
            HLA <- prepareMiDAS(
            hla_calls = dat_HLA,
            colData = pheno,
            experiment = "hla_alleles"
            )

            HLA_model <- glm(reformulate(covariates, response = phecode), 
                            data = HLA, family = binomial())
            HLA_results <- runMiDAS(
            object = HLA_model, 
            experiment = "hla_alleles", 
            inheritance_model = "dominant",
            lower_frequency_cutoff = 0.01, 
            upper_frequency_cutoff = 0.99, 
            correction = "bonferroni", 
            exponentiate = TRUE
            )
            samplenum <- dim(dat_HLA)[1]
            filename = paste(phecode, ".", samplenum, ".hla_alleles.csv", sep = "")
            write.csv(HLA_results, file = filename, row.names = FALSE)

            # hla - aa association
            # AA fine mapping
            HLA_AA <- tryCatch({
              prepareMiDAS(
                hla_calls = dat_HLA,
                colData = pheno,
                experiment = "hla_aa"
              )
            }, error = function(e) {
              cat("Error creating HLA_AA object:", e$message, "\n")
              return(NULL)
            })

            # Only proceed if HLA_AA was successfully created
            if (!is.null(HLA_AA) && (inherits(HLA_AA, "SummarizedExperiment") || is.data.frame(HLA_AA) || is.matrix(HLA_AA))) {
              cat("HLA_AA object successfully created. Proceeding with analysis.\n")
              HLA_AA_model <- glm(reformulate(covariates, response = phecode), 
                                  data = HLA_AA, family = binomial())
              HLA_AA_omnibus_results <- runMiDAS(
                HLA_AA_model,
                experiment = "hla_aa",
                inheritance_model = "dominant",
                conditional = FALSE,
                omnibus = TRUE,
                lower_frequency_cutoff = 0.01,
                upper_frequency_cutoff = 0.99,
                correction = "bonferroni"
              )
              aa_filename = paste(phecode, ".", samplenum, ".hla_aa.csv", sep = "")
              write.csv(HLA_AA_omnibus_results, file = aa_filename, row.names = FALSE)
            } else {
              cat("Warning: HLA_AA object creation failed or returned invalid object. Skipping HLA-AA analysis.\n")
              aa_filename = paste(phecode, ".", samplenum, ".hla_aa.csv", sep = "")
              write.csv(data.frame(message = "HLA-AA analysis skipped - prepareMiDAS failed"), 
                        file = aa_filename, row.names = FALSE)
            }


            # KIR gene association
            midas <- prepareMiDAS(
            hla_calls = dat_HLA,
            kir_calls = kir_calls,
            colData = pheno,
            experiment = c("hla_NK_ligands","kir_genes", "hla_kir_interactions")
            )
            KIR_model <- glm(reformulate(covariates, response = phecode), 
                            data = midas, family = binomial())
            KIR_results <- runMiDAS(
            KIR_model,
            experiment = "kir_genes",
            lower_frequency_cutoff = 0.01,
            upper_frequency_cutoff = 0.99,
            exponentiate = TRUE
            )
            kir_filename = paste(phecode, ".", samplenum, ".kir_gene.csv", sep = "")
            write.csv(KIR_results, file = kir_filename, row.names = FALSE)

            HLA_NK_results <- runMiDAS(
            KIR_model,
            experiment = "hla_NK_ligands",
            inheritance_model = "dominant",
            lower_frequency_cutoff = 0.01,
            upper_frequency_cutoff = 0.99,
            exponentiate = TRUE
            )
            hla_nk_filename = paste(phecode, ".", samplenum, ".hla_nk.csv", sep = "")
            write.csv(HLA_NK_results, file = hla_nk_filename, row.names = FALSE)

            HLA_KIR_results <- runMiDAS(
            KIR_model,
            experiment = "hla_kir_interactions",
            lower_frequency_cutoff = 0.01,
            upper_frequency_cutoff = 0.99,
            exponentiate = TRUE
            )
            hla_kir_filename = paste(phecode, ".", samplenum, ".hla_kir.csv", sep = "")
            write.csv(HLA_KIR_results, file = hla_kir_filename, row.names = FALSE)

        # Plot
            HLA_results["gene"] = sapply(HLA_results$allele, function(a) strsplit(a, "*", fixed = TRUE)[[1]][1])
            
            hla_plot_obj <- ggplot(HLA_results[1:10,], aes(x = estimate, y = reorder(allele, -p.value), color = gene)) +
            geom_point(size = 3) +
            xlim(0, 20) + 
            geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, orientation = "y") +
            geom_vline(xintercept = 1, linetype = "dashed") + # The "No Effect" line
            theme_minimal() + 
            labs(title = phecode, x = "Odds Ratio", y = "Alleles")
            hla_plot = paste(phecode, "_", samplenum, "_hla_allele.png", sep = "")
            ggsave(hla_plot, plot = hla_plot_obj, width = 6, height = 4, units = "in")

            #Plot Manhattan
            sorted_HLA_results <- HLA_results %>% arrange(allele)
            sorted_HLA_results["gene"] = sapply(sorted_HLA_results$allele, function(a) strsplit(a, "*", fixed = TRUE)[[1]][1])
            sorted_HLA_results <- sorted_HLA_results %>% group_by(gene) %>% mutate(position = row_number()) %>% ungroup() # Optional: ungroup the data after the operation

            plot_data <- sorted_HLA_results %>%
            mutate(
                # Calculate -log10 P-value
                logP = -log10(p.value),
                # Create a continuous X-axis for multiple genes
                # (This ensures HLA-B plots to the right of HLA-A, not on top of it)
                bp_cum = position + 
                case_when(gene == "A" ~ 0,
                            gene == "B" ~ 50,
                            gene == "C" ~ 100,
                            gene == "DPA1" ~ 150,
                            gene == "DPA2" ~ 200,
                            gene == "DPB1" ~ 250,
                            gene == "DPB2" ~ 300,
                            gene == "DQA1" ~ 350,
                            gene == "DQA2" ~ 400,
                            gene == "DQB1" ~ 450,
                            gene == "DQB2" ~ 500,
                            gene == "DRA" ~ 550,
                            gene == "DRB1" ~ 600,
                            gene == "DRB3" ~ 650,
                            gene == "DRB4" ~ 700,
                            gene == "DRB5" ~ 750, 
                            gene == "DPB2" ~ 800, 
                            gene == "DQB1" ~ 850),
                label = ifelse(logP > 2, as.character(allele), "")
            )
            
            # Calculate center positions for x-axis labels
            axis_set <- plot_data %>% 
                group_by(gene) %>% 
                summarize(center = mean(bp_cum))

            # 3. Generate the Plot
            manhattan_plot_obj <- ggplot(plot_data, aes(x = bp_cum, y = logP, color = gene)) +
            # The points
            geom_point(alpha = 0.8, size = 1.5) +
            # The significance threshold line (e.g., Bonferroni)
            geom_hline(yintercept = 2, color = "red", linetype = "dashed") +
            # Formatting axes
            scale_x_continuous(label = axis_set$gene, breaks = axis_set$center) +
            # add labels
            geom_text_repel(
                aes(label = label),
                box.padding = 0.5, # Adjust padding around the label
                point.padding = 0.5 # Adjust padding around the point
            ) + 
            # Labels and Theme
            labs(title = phecode,
                x = "HLA Alleles",
                y = "-log10(P-value)") +
            theme_classic() +
            theme(
                legend.position = "none",                # Hide legend if labels are on X-axis
                axis.text.x = element_text(size = 8, face = "bold"),
                panel.grid.major.y = element_line(color = "grey90")
            )
            manhattan_plot = paste(phecode, "_", samplenum, "_hla_manhattan.png", sep = "")
            ggsave(manhattan_plot, plot = manhattan_plot_obj, width = 6, height = 4, units = "in")
            ' > run_midasHLA.R

    Rscript run_midasHLA.R ~{hla_alleles_types} ~{kir_alleles_types} ~{phenotype_file} ~{phecode} ~{hla_resolution} > midashla.log 2>&1
    >>>

  output {
    Array[File] outputfiles = glob("~{phecode}*")
    Array[File] outputfigures = glob("*.png")
    File midasHLA_log = "midashla.log"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/midashla:v1"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }

}