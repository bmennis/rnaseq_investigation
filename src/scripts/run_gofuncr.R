#Load gofuncr for analysis
library(GOfuncR)
library(data.table)
library(tidyverse)

#load deseq2 output for filtering and analysis, set as tibble, and group by symbol
deseq2_out <- read.csv("Desktop/deseq2_output_10_12_18/condition_treated_results_filtered_tx2gene_full_mouse.csv")
deseq2_out <- as.tibble(deseq2_out)
deseq2_out <- deseq2_out %>% group_by(symbol)

#filter out NA from output columns
deseq2_out_filtered <- deseq2_out[complete.cases(deseq2_out[,5:8]),]

#filter output to remove duplicate gene symbols selecting for lowest padj value and one row of data for exact matches
deseq2_out_filtered <- deseq2_out_filtered %>% 
                                          filter(padj == min(padj)) %>%
                                          filter(row_number() == 1)

#take the filtered deseq results with padj < 0.001 for candidate and set all else as background
candidate_genes <- deseq2_out_filtered %>% filter(padj < 0.001)
background_genes <- deseq2_out_filtered %>% filter(padj >= 0.001)

#make a data frame of cand and bg genes by selecting appropriate symbols
mice_cand_genes = data.frame(candidate_genes$symbol)
mice_bg_genes = data.frame(background_genes$symbol)

#convert the data frames to matrices to make one long list of cand followed by bg
mice_bg_genes <- as.matrix(mice_bg_genes)
mice_cand_genes <- as.matrix(mice_cand_genes)
mice_genes <- c(mice_cand_genes,mice_bg_genes)

#make a list of 1 and 0 for cand and bg respectively for each gene symbol list length to match with symbols
mice_gene_types <- c(rep(1,length(mice_cand_genes)),rep(0,length(mice_bg_genes)))

#make data frame of one column cand and bg gene symbols and another of 1 for cand and 0 for bg
mice_gene_dat <- data.frame(mice_genes,mice_gene_types)

#run gofuncr on data frame of cand and bg genes with mouse data base list to get enrichment
res_hyper_mouse = go_enrich(mice_gene_dat, organismDb='Mus.musculus')

#write the 3 results of the test to csv files
write_csv(data.frame(res_hyper_mouse[1]), "Desktop/gofuncr_results.csv")
write_csv(data.frame(res_hyper_mouse[2]), "Desktop/gofuncr_results_input.csv")
write_csv(data.frame(res_hyper_mouse[3]), "Desktop/gofuncr_results_anno.csv")

#write deseq2 filtered output to csv file for GSEA analysis
write_csv(deseq2_out_filtered, "Desktop/deseq2_results_filtered.csv")

