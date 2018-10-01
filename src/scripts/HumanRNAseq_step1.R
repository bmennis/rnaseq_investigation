
##### pre-1) change the workinng directory where kallisto results stored
setwd("/mnt/isilon/dbhi_bfx/xie/Brian_folder/rnaseq_investigation/data/raw/fq_files/")

dir<-"/mnt/isilon/dbhi_bfx/xie/Brian_folder/rnaseq_investigation/data/raw/fq_files/"
##### pre-2) load deseq2
library(DESeq2)
library(dplyr)

library("tximport")
library("readr")
library(rhdf5)

samples <- read.table("/mnt/isilon/dbhi_bfx/xie/Brian_folder/rnaseq_investigation/data/interim/annotations/directory_sample.txt", header = TRUE)

files <- file.path(dir, samples$sample, "abundance.h5")
names(files)<-samples$sample
tx2gene <- read.table("/mnt/isilon/dbhi_bfx/xie/Brian_folder/rnaseq_investigation/data/interim/annotations/tx2gene_mouse.txt", header = T)
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene,ignoreTxVersion=T, txOut = F)
head(txi.kallisto$counts)

dds <- DESeqDataSetFromTximport(txi.kallisto,
                                colData = samples,
                                design = ~ condition)

keep <- rowSums(counts(dds)) >0
dds <- dds[keep,]
dds$condition <- factor(dds$condition)

dds <- DESeq(dds)
