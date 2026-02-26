version 1.0

workflow PheWAS{
    meta{
        description: "a workflow for haplotype PheWAS association"
    }
    input{
        File? genotype_file
        File phenotype_file
        File meta_data_file
        File? genotype_json
        File phecodex_description
        Int allele_resolution
        Boolean additive_mode

        String memory = "8G"
        Int cpu = 2
        Int disk_size = 50
        Int preemptible = 1

    }

    if (defined(genotype_json)) {
        call PreProcessFileFromJSON{input:
            meta_data_file = meta_data_file,
            phecodex_info = phenotype_file,
            genotype_json = select_first([genotype_json,genotype_file]),
            allele_resolution = allele_resolution,
            disk_size = disk_size
        }
    }
    if (defined(genotype_file)) {
        call PreProcessFile{input:
            meta_data_file = meta_data_file,
            phecodex_info = phenotype_file,
            genotype_info = select_first([genotype_json,genotype_file]),
            disk_size = disk_size
        }
    }

    scatter (haplotype in select_first([PreProcessFileFromJSON.haplotype_list,PreProcessFile.haplotype_list])) {
        call RunPheWAS {input:
            genotype_file =select_first([PreProcessFileFromJSON.genotype_file,PreProcessFile.genotype_file]),
            phenotype_file =select_first([PreProcessFileFromJSON.phenotype_file,PreProcessFile.phenotype_file]),
            meta_data_file = select_first([PreProcessFileFromJSON.metadata_file,PreProcessFile.metadata_file]),
            haplotype = haplotype,
            additive_genotype = additive_mode,
            memory = memory,
            cpu = cpu,
            disk_size = disk_size,
            preemptible = preemptible
        }        
    }
    
    output {
        Array[File] PheWas_results = flatten(RunPheWAS.outputfiles)
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

task PreProcessFileFromJSON {
    input {
        File meta_data_file
        File phecodex_info
        File genotype_json
        Int allele_resolution
        Int disk_size
        Int minimal_case_num = 20
        Float minimal_af = 0.0001

        RuntimeAttr? runtime_attr_override

    }

    command <<<
        set -euxo pipefail

        python - --phecodex ~{phecodex_info} \
                --genotype ~{genotype_json} \
                --metadata ~{meta_data_file} \
                --case_num ~{minimal_case_num} \
                --threshold ~{minimal_af} \
                --resolution ~{allele_resolution} \
                <<-'EOF'
        import os
        import pandas as pd
        import subprocess
        import numpy
        from datetime import date
        import argparse
        import json


        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--phecodex',
                                type=str)
            parser.add_argument('--genotype',
                                type=str)
            parser.add_argument('--metadata',
                                type=str)
            parser.add_argument('--case_num',
                    type=int)
            parser.add_argument('--threshold',
                    type=float)
            parser.add_argument('--resolution',
                    type=int)

            args = parser.parse_args()

            # load metadata info
            df_meta = pd.read_csv(args.metadata, sep = ",", header=0, dtype=str)
            print(df_meta.head())
            
            # load genotypes

            threshold = args.threshold
            with open(args.genotype, "r") as fp:
                D = json.load(fp)

            Info = {}
            resolution = args.resolution
            pid_set = set()
            for haplotype, l in D.items():
                if haplotype.startswith("KIR"):
                    if resolution == 0:
                        genelist = []
                        for h in haplotype.split(","):
                            genelist.append(h.split("*")[0])
                    elif resolution == 1:
                        genelist = []
                        for h in haplotype.split(","):
                            genebody = h.split("*")[0]
                            allele = h.split("*")[1][:3]
                            genelist.append(genebody + "*" + allele)                  
                    elif resolution == 2:
                        genelist = []
                        for h in haplotype.split(","):
                            genebody = h.split("*")[0]
                            allele = h.split("*")[1][:5]
                            genelist.append(genebody + "*" + allele)   
                    else:
                        genelist = haplotype.split(",")              
                else:
                    genelist = [":".join(h.split(":")[:resolution]) for h in haplotype.split(",")]
                hap = ",".join(genelist)
                Info[hap] = Info.get(hap, {})
                for hapid in l:
                    pid = hapid.split("_")[0]
                    pid_set.add(pid)
                    Info[hap][pid] = Info[hap].get(pid, 0) + 1
            df_genotype_all = pd.DataFrame(Info).fillna(0)
            df_genotype_all.index = df_genotype_all.index.astype(str)

            df_freq = df_genotype_all.sum(axis = 0) / (df_genotype_all.shape[0] * 2)
            df_genotype = df_genotype_all.loc[:, df_freq > threshold]
            df_genotype.index = df_genotype.index.astype(str)
            print(df_genotype.head())

            # load phecodex
            filter_pheno_df = pd.read_csv(args.phecodex, header=0, index_col=0, dtype=str)
            # replace TRUE, FALSE with 1 and 0
            replacement_map = {'True': 1, 'False': 0, "FALSE": 0, "TRUE": 1, "NaN": numpy.nan}

            # Replace the strings with numbers
            filter_pheno_df = filter_pheno_df.replace(replacement_map)

            # metadata filtering
            df_meta = df_meta[['person_id', 'age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'sex_at_birth']].drop_duplicates()
            indexlist = set(filter_pheno_df['person_id'].astype(str).tolist()) & set(df_meta['person_id'].astype(str).tolist()) & set(df_genotype.index.astype(str).tolist())

            print(len(indexlist))
            indexlist = sorted(indexlist)
            df_meta = df_meta.loc[df_meta['person_id'].isin(indexlist), :]
            df_meta = df_meta.set_index("person_id")
            filter_pheno_df = filter_pheno_df.loc[filter_pheno_df['person_id'].isin(indexlist), :]
            filter_pheno_df = filter_pheno_df.set_index("person_id")
            df_genotype = df_genotype.loc[indexlist, :]
            df_genotype['person_id'] = df_genotype.index.tolist()

            # filter phecodex data given selected cohort
            col_sums = filter_pheno_df.astype(float).sum(axis = 0)
            condition = col_sums > args.case_num
            filter_pheno_df = filter_pheno_df.loc[:,condition]

            df_pheno_new = filter_pheno_df.astype(float) > 0.5

            replacement_map_sex = {'PMI: Skip': numpy.nan, 
                   'Sex At Birth: Sex At Birth None Of These': numpy.nan, 
                   'I prefer not to answer':numpy.nan,
                   'Intersex' : numpy.nan
                  }
            df_meta = df_meta.replace(replacement_map_sex)

            df_genotype.astype(float).to_csv("genotypes.csv", index = True)
            df_pheno_new.to_csv("phenotypes.csv", index=True)
            df_meta.to_csv("metadata.csv", index=True)

            with open("haplotypelist.txt", 'w') as fp:
                for item in df_genotype.columns.tolist():
                    if item == "person_id":
                        continue
                    fp.writelines(item + "\n")

        if __name__ == "__main__":
            main()

        EOF

    >>>

    output {
        File genotype_file = "genotypes.csv"
        File phenotype_file = "phenotypes.csv"
        File metadata_file = "metadata.csv"
        Array[String] haplotype_list = read_lines("haplotypelist.txt")
        
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/midashla:v0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task PreProcessFile {
    input {
        File meta_data_file
        File phecodex_info
        File genotype_info
        Int disk_size
        Int minimal_case_num = 20

        RuntimeAttr? runtime_attr_override

    }

    command <<<
        set -euxo pipefail

        python - --phecodex ~{phecodex_info} \
                --genotype ~{genotype_info} \
                --metadata ~{meta_data_file} \
                --case_num ~{minimal_case_num} \
                <<-'EOF'
        import os
        import pandas as pd
        import subprocess
        import numpy
        from datetime import date
        import argparse


        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--phecodex',
                                type=str)
            parser.add_argument('--genotype',
                                type=str)
            parser.add_argument('--metadata',
                                type=str)
            parser.add_argument('--case_num',
                    type=int)

            args = parser.parse_args()

            # load metadata info
            df_meta = pd.read_csv(args.metadata, sep = ",", header=0, dtype=str)
            print(df_meta.head())
            
            # load genotypes
            df_genotype = pd.read_csv(args.genotype, index_col = 0, dtype=str)
            df_genotype.index = df_genotype.index.astype(str)
            print(df_genotype.head())

            # load phecodex
            filter_pheno_df = pd.read_csv(args.phecodex, header=0, index_col=0, dtype=str)
            # replace TRUE, FALSE with 1 and 0
            replacement_map = {'True': 1, 'False': 0, "FALSE": 0, "TRUE": 1, "NaN": numpy.nan}

            # Replace the strings with numbers
            filter_pheno_df = filter_pheno_df.replace(replacement_map)

            # metadata filtering
            df_meta = df_meta[['person_id', 'age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'sex_at_birth']].drop_duplicates()
            indexlist = set(filter_pheno_df['person_id'].astype(str).tolist()) & set(df_meta['person_id'].astype(str).tolist()) & set(df_genotype.index.astype(str).tolist())

            print(len(indexlist))
            indexlist = sorted(indexlist)
            df_meta = df_meta.loc[df_meta['person_id'].isin(indexlist), :]
            df_meta = df_meta.set_index("person_id")
            filter_pheno_df = filter_pheno_df.loc[filter_pheno_df['person_id'].isin(indexlist), :]
            filter_pheno_df = filter_pheno_df.set_index("person_id")
            df_genotype = df_genotype.loc[indexlist, :]
            df_genotype['person_id'] = df_genotype.index.tolist()

            # filter phecodex data given selected cohort
            col_sums = filter_pheno_df.astype(float).sum(axis = 0)
            condition = col_sums > args.case_num
            filter_pheno_df = filter_pheno_df.loc[:,condition]

            df_pheno_new = filter_pheno_df.astype(float) > 0.5

            replacement_map_sex = {'PMI: Skip': numpy.nan, 
                   'Sex At Birth: Sex At Birth None Of These': numpy.nan, 
                   'I prefer not to answer':numpy.nan,
                   'Intersex' : numpy.nan
                  }
            df_meta = df_meta.replace(replacement_map_sex)

            df_genotype.astype(float).to_csv("genotypes.csv", index = True)
            df_pheno_new.to_csv("phenotypes.csv", index=True)
            df_meta.to_csv("metadata.csv", index=True)

            with open("haplotypelist.txt", 'w') as fp:
                for item in df_genotype.columns.tolist():
                    if item == "person_id":
                        continue
                    fp.writelines(item + "\n")

        if __name__ == "__main__":
            main()

        EOF

    >>>

    output {
        File genotype_file = "genotypes.csv"
        File phenotype_file = "phenotypes.csv"
        File metadata_file = "metadata.csv"
        Array[String] haplotype_list = read_lines("haplotypelist.txt")
        
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/midashla:v0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task RunPheWAS {

  input {
    File genotype_file
    File phenotype_file
    File meta_data_file
    String haplotype
    Boolean additive_genotype

    # Runtime parameters
    String memory
    Int cpu
    Int disk_size
    Int preemptible

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euxo pipefail

    echo '
            args <- commandArgs(trailingOnly = TRUE)
            genotype_file <- args[1]
            phenotype_file <- args[2]
            meta_data_file <- args[3]
            haplotype <- args[4]
            additive_mode <- args[5]
            
            library(tidyverse)
            library(dplyr)
            library(parallel)
            library(PheWAS)
            library(readr)

            # load genotype
            genotypes <- read_csv(genotype_file)
            colnames(genotypes)[1] <- "person_id"
            print(genotypes)

            # load phenotype
            phenotypes <- read_csv(phenotype_file)

            # load metadata
            metadata <- read_csv(meta_data_file)

            # perform PheWAS

            print(haplotype)
            cols_to_select <- c("person_id", haplotype)
            selected_cols <- genotypes[, cols_to_select]
            colnames(selected_cols)[1] <- "person_id"

            data <- inner_join(phenotypes, selected_cols, by = "person_id")
            data <- inner_join(data, metadata, by = "person_id")
            print(data)
            results <- phewas_ext(data,
                            phenotypes=names(phenotypes)[c(-1)],  # All phecode columns
                            genotypes=c(haplotype),
                            covariates=c("sex_at_birth","age","PC1",
                                        "PC2","PC3", "PC4", "PC5"), 
                            additive.genotypes = additive_mode, 
                            cores=8)
            
            results <- results[!is.na(results$p), ]
            results$p_bonferroni <- p.adjust(results$p, method = "bonferroni")

            #List the significant results
            sig_results <- results[results$p_bonferroni<0.05,]
            num_of_sig = dim(sig_results)[1]
            print(c(haplotype,num_of_sig))


            haplotype_clean <- gsub("KIR", "", haplotype)
            haplotype_clean_ <- gsub("HLA", "", haplotype_clean)
            output_file <- sprintf("phewas_significant_%s_results.csv", num_of_sig)
            write.csv(results, file = output_file, row.names = FALSE) 
            
                     
        ' > run_PheWAS.R

    Rscript run_PheWAS.R ~{genotype_file} ~{phenotype_file} ~{meta_data_file} ~{haplotype} ~{additive_genotype} > phewaslog.log 2>&1
    >>>

  output {
    Array[File] outputfiles = glob("*_results.csv")
    File logfile = "phewaslog.log"
  }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/phewas:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

}


task PlotResults {
    input {
        File phewas_results
        File phecodex_info
        Int disk_size

        RuntimeAttr? runtime_attr_override

    }

    command <<<
        set -euxo pipefail

        pip install adjustText

        python - --phewas_results ~{phewas_results} \
                --phecodex_info ~{phecodex_info} \
                <<-'EOF'
        import os
        import argparse
        import pandas as pd
        import subprocess
        import numpy
        from adjustText import adjust_text
        import matplotlib.pyplot as plt
        from matplotlib import colormaps

        def plot_phewas_manhattan(filename, category_description,phecode_discription, annotation = True):
            df_results = pd.read_csv(filename, index_col = False)
            
            if (df_results['p'] == 0).all() or (df_results['p'].isna()).all():
                print("Warning", filename)
                return

            df_results['logP'] = [round(-numpy.log10(p), 2) if p > 0 else 10 for p in df_results['p'].tolist()]

            df_results["phecode_category"] = [category_description[p.split("_")[0]] for p in df_results['phenotype'].tolist()]

            df_results.phecode_category = df_results.phecode_category.astype('category')
            #df_results.phecode_category = df_results.phecode_category.cat.set_categories(['ch-%i' % i for i in range(12)], ordered=True)
            df_results = df_results.sort_values('phecode_category')

            # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
            df_results['ind'] = range(len(df_results))
            df_grouped = df_results.groupby(('phecode_category'))

            # manhattan plot
            fig = plt.figure(figsize=(14, 8)) # Set the figure size
            ax = fig.add_subplot(111)
            #colors = ['pink','lightgreen','lightblue', 'gold']
            colors = colormaps['Set2'].colors
            x_labels = []
            x_labels_pos = []
            ypos = set()
            texts = []
            for num, (name, group) in enumerate(df_grouped):
                group.plot(kind='scatter', x='ind', y='logP',color=colors[num % len(colors)], ax=ax)
                x_labels.append(name)
                x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
                if annotation:
                    df_sub = group[group['p_bonferroni']<0.05]
                    #print(df_sub)
                    for i, row in df_sub.iterrows():
                        phecode = row['phenotype']
                        txt = phecode + ":" + phecode_discription[phecode]
                        x = row['ind']
                        y = row['logP']
                        texts.append(ax.text(x, y, txt, fontsize=10))
                        
                        
            ax.set_xticks(x_labels_pos)
            ax.set_xticklabels(x_labels, rotation = 90)
            # add text
            #adjust_text(texts, arrowprops=dict(arrowstyle="->", color='black', lw=1))
            adjust_text(texts)

            # set axis limits
            ax.set_xlim([0, len(df_results)])
            ax.set_ylim([0, df_results['logP'].max()+1])

            # x axis label
            ax.set_xlabel('PhecodeX', fontsize = 15)
            ax.set_ylabel('LogP', fontsize = 15)
            
            # Get title from filename
            title = filename.split("/")[-1].split("_")[0]
            plt.title(title, fontsize=10, fontweight='bold')
            
            # save figure
            output_filename = filename.replace(".csv", ".pdf")
            plt.tight_layout()
            fig.savefig(output_filename, dpi=300, bbox_inches='tight', format='pdf')
            print(f"Figure saved to: {output_filename}")
            plt.show() 
            
            return ax
    


        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--phewas_results',
                                type=str)
            parser.add_argument('--phecodex_info',
                                type=str)

            args = parser.parse_args()
            results_file = args.phewas_results
            phecodex_file = args.phecodex_info

            # load Phecodex info
            phe_index = pd.read_csv(phecodex_file)
            phecode_discription = {}
            category_description = {}
            for index, row in phe_index.iterrows():
                phecode = row['phecodex']
                description = row['description']
                group = row['group']
                
                category = phecode.split("_")[0]
                phecode_discription[phecode] = description
                category_description[category] = group

            # plot PheWas
            ax = plot_phewas_manhattan(results_file, category_description, phecode_discription, annotation = True)

        if __name__ == "__main__":
            main()

        EOF

    >>>

    output {
        Array[File]? outputfiles = glob("*_results.pdf")
        
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/midashla:v0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}