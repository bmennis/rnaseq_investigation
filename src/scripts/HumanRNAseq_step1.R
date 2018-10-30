args = commandArgs(trailingOnly=TRUE)

in_dir <- args[1]
sample_file <- args[2]
gene_input <- args[3]
out_csv <- args[4]
out_pca <- args[5]
out_cluster <- args[6]

##### pre-1) change the workinng directory where kallisto results stored
#setwd("/mnt/isilon/dbhi_bfx/xie/Brian_folder/rnaseq_investigation/data/interim/deseq2_output/")

#dir<-"/mnt/isilon/dbhi_bfx/xie/Brian_folder/rnaseq_investigation/data/raw/fq_files/"
##### pre-2) load deseq2
library(DESeq2)
library(dplyr)
library("AnnotationDbi")
library("org.Mm.eg.db")
library("tximport")
library("readr")
library(rhdf5)
library("vsn")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("IHW")
library("ggrepel")

samples <- read.table(sample_file, header = TRUE)

files <- file.path(in_dir, samples$sample, "abundance.h5")
names(files)<-samples$sample
tx2gene <- read.table(gene_input, header = T)
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene,ignoreTxVersion=T, txOut = F)
head(txi.kallisto$counts)

dds <- DESeqDataSetFromTximport(txi.kallisto,
                                colData = samples,
                                design = ~ condition)

keep <- rowSums(counts(dds)) >0
dds <- dds[keep,]
dds$condition <- factor(dds$condition)

dds <- DESeq(dds)

res <- results(dds)

res$symbol <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")


resOrdered <- res[order(res$pvalue),]

write.csv(as.data.frame(resOrdered), 
          file= out_csv)

colData(dds)

vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
head(assay(vsd), 3)

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sample, sep=":")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(out_cluster)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

pdf(out_pca)
pcaData <- plotPCA(vsd, intgroup=c("condition", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaPlot <- ggplot(pcaData, aes(PC1, PC2, color=condition, label="")) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  guides(label = TRUE)
pcaPlot + 
  geom_label_repel(aes(label = sample),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50')
dev.off()
