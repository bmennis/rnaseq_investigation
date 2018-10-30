#load deseq2 data to filter and organize to run in gsea
deseq2_filtered <- read.csv("Desktop/deseq2_results_filtered.csv")

#load annotated mouse gene symbols with homologous human 
#gene symbols
human_mouse_homologs <- read.table("Desktop/human_mouse_homolog.out", header = TRUE)
human_mouse_homologs <- as.data.frame(human_mouse_homologs)

#create data frame of gene symbol column and stat column 
#from deseq2
gsea_dat <- data.frame(deseq2_filtered[,"symbol"],deseq2_filtered[,"stat"])
colnames(gsea_dat) <- c('id','stat')

#filter data with human_mouse_homologs to match stat 
#with corresponding homolog
#start by merging gsea_dat with the homologs by gene 
#symbols to get HomoloGene_ID
gsea_homologene_ids <- merge(gsea_dat, human_mouse_homologs, by.x = "id", by.y = "Symbol")[,1:3]

#next merge the homologene id data back with the homolog data
#by homologene id to get both mouse and human symbols but only
#for the data from the input
gsea_mouse_human_dat <- merge(gsea_homologene_ids, human_mouse_homologs, by.x = "HomoloGene_ID", by.y = "HomoloGene_ID")

#last filter out all rows for mouse so only human data is left
#and take just the 'Symbol' and 'stat' columns for the rnk file
gsea_human_dat <- gsea_mouse_human_dat %>% filter(Common_Organism_Name != "mouse")
mice_human_homolog_gsea_input <- data.frame(Symbol = gsea_human_dat[,'Symbol'],Stat = gsea_human_dat[,'stat'])

#write data frame to tsv file named .rnk so that it can run in gsea 
write_tsv(mice_human_homolog_gsea_input,"Desktop/mice_human_homolog_gsea_input.rnk")














