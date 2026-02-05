version 1.0

workflow PheWAS{
    meta{
        description: "a workflow for haplotype PheWAS association"
    }
    input{
        File genotype_file
        File phenotype_file
        File meta_data_file
        File phecodex_description
        String haplotype

        String memory = "8G"
        Int cpu = 2
        String disk_size = "local-disk 50 HDD"
        Int preemptible = 1

    }


    call RunPheWAS {input:
        genotype_file = genotype_file,
        phenotype_file = phenotype_file,
        meta_data_file = meta_data_file,
        haplotype = haplotype,
        memory = memory,
        cpu = cpu,
        disk = disk_size,
        preemptible = preemptible
    }

    call PlotResults{input:
        phewas_results = select_first(RunPheWAS.outputfiles),
        phecodex_info = phecodex_description,
        disk_size = disk_size
    }

    output {
        Array[File] PheWas_results = RunPheWAS.outputfiles
        Array[File] figures = PlotResults.outputfiles
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

task RunPheWAS {

  input {
    File genotype_file
    File phenotype_file
    File meta_data_file
    String haplotype

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
            genotype_file <- args[1]
            phenotype_file <- args[2]
            meta_data_file <- args[3]
            haplotype <- args[4]


            install.packages("devtools")
            install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
            devtools::install_github("PheWAS/PheWAS")
            
            library(tidyverse)
            library(dplyr)
            library(parallel)
            library(PheWAS)
            library(readr)

            # load genotype
            genotypes <- read_csv(genotype_file)
            colnames(genotypes)[1] <- "person_id"

            # load phenotype
            phenotypes <- read_csv(phenotype_file)

            # load metadata
            metadata <- read_csv(meta_data_file)

            # perform PheWAS

            print(haplotype)
            cols_to_select <- c("person_id", haplotype)
            selected_cols <- genotypes[, cols_to_select]
            colnames(selected_cols)[1] <- 'person_id'
            data <- inner_join(phenotypes, selected_cols, by = "person_id")
            data <- inner_join(data, metadata, by = "person_id")
            results <- phewas_ext(data,
                            phenotypes=names(phenotypes)[-1],  # All phecode columns
                            genotypes=c(haplotype),
                            covariates=c("sex_at_birth",'age','PC1',
                                        'PC2','PC3', 'PC4', "PC5"), 
                            cores=8)
            
            n_tests <- nrow(results[!is.na(results$p), ])
            results$p_bonferroni <- p.adjust(results$p, method = "bonferroni")

            #List the significant results
            sig_results <- results[results$p_bonferroni<0.05,]
            num_of_sig = dim(sig_results)[1]
            print(c(haplotype,num_of_sig))
            if (num_of_sig > 0){
                output_file <- sprintf('%s_significant_%s_results.csv', haplotype, num_of_sig)
                write.csv(results, file = output_file, row.names = FALSE) 
            }
                
        ' > run_PheWAS.R

    Rscript run_PheWAS.R ~{genotype_file} ~{phenotype_file} ~{meta_data_file} ~{haplotype} > phewaslog.log 2>&1
    >>>

  output {
    Array[File] outputfiles = glob("*_results.csv")
    File logfile = "phewaslog.log"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/midashla:v1"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }

}

task PlotResults {
    input {
        File phewas_results
        File phecodex_info
        Int disk_size

    }

    command <<<
        set -euxo pipefail

        pip install adjustText

        python - --phewas_results ~{phewas_results}                 --phecodex_info ~{phecodex_info}                 <<-'EOF'
        import os
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
        Array[File] outputfiles = glob("*_results.pdf")
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/midashla:v0"
        memory: "4 GB"
        cpu: 1
        disks: disk_size
    }
}
