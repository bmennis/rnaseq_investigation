#######################################
#Code to merge tsv file and generate FPKM file
#######################################

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))

inDir <- args[1]
outFile <- args[2]
rm(args)

#Read in all TSV Files
setwd(inDir);
tsvFiles <- list.files()[grep("tsv", list.files())];

#Remove files with no values
isEmpty <- function(x)
{
cha <- file.info(x)[[1]]>0
return(cha);
}
tsvFiles <- tsvFiles[as.vector(sapply(tsvFiles, FUN=isEmpty))];

tpmDF <- read.delim(tsvFiles[1], header=T);
tpmDF <- tpmDF[,c(1,5)];
colnames(tpmDF)<- c("Transcript", gsub(".tsv", "", tsvFiles[1]));

for(i in 2:length(tsvFiles))
{
tmp <- read.delim(tsvFiles[i], header=T);
tpmDF <- cbind(tpmDF, tmp[5]);
colnames(tpmDF)[i+1]<- gsub(".tsv", "", tsvFiles[i]);
}

rownames(tpmDF) <- tpmDF[,1];
tpmDF <- tpmDF[-1];

write.table(tpmDF, outFile, sep="\t", row.names=T);
