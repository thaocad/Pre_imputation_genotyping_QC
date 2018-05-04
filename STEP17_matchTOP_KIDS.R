library(readr)
library(dplyr)
#setwd("//cortex.meaney.lab/tnguyen/NFP_QC/STEP2-QC/STEP17_TOP")
args = commandArgs(trailingOnly=TRUE)

#kidsinputmap=read_delim("STEP17-TOP/LEARS_09Mar2018_step17.map", col_names = F,delim = "\t" )
kidsinputmap=read_delim(args[1], col_names = F,delim = "\t" )
#kidsinputgen=read_delim("STEP17-TOP/LEARS_09Mar2018_step17.gen", col_names = F,delim = " " )
kidsinputgen=read_delim(args[2], col_names = F,delim = " " )

#topArrayLines=read_lines("RESOURCES/GSA-24v1-0_A2.update_alleles.txt")
topArrayLines=read_lines(args[3])
topArrayLinesSub <- gsub(pattern=" ",replacement="\t",topArrayLines,fixed=TRUE)
#strsplit, make dataframe
topArray <- do.call(rbind.data.frame,strsplit(topArrayLinesSub,split="\t"))

colnames(kidsinputmap)=c("CHROM", "SNP", "MORGANS", "BP")
colnames(kidsinputgen)[1:5]=c("CHROM", "SNP", "BP", "A1", "A2")
colnames(topArray)=c("SNP", "A","B", "TOP1", "TOP2")
topArray=apply(topArray, 2, as.character)

## see if the snps have the same allels as in the TOP file
kidsInfo=kidsinputgen[1:5]
kidsInfoTop=inner_join(kidsInfo, as.data.frame(topArray), by="SNP" , suffix=c("KIDS", "TOP") , copy=T)

### see if the non ATCG file 
topArrayA1=topArray[,c("SNP","TOP1")]
write.table(topArrayA1, "STEP17-TOP/Array_TOP1.txt", quote = F, row.names = F, col.names = F)

topArrayA2=topArray[,c("SNP","TOP2")]
write.table(topArrayA2, "STEP17-TOP/Array_TOP2.txt", quote = F, row.names = F, col.names = F)
