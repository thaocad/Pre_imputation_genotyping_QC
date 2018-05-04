library(data.table)
library(tidyverse)
library(magrittr)
args = commandArgs(trailingOnly=TRUE)
# requirement : GEN file aligned in TOP-strand
# path = complete path of GEN file without the extension, i.e., sans ".gen", outputs will also be recorded here, unless specified with outpath

remove_duptrip_indels_23to26 <- function(filepath, outpath=NULL){
  if(is.null(outpath)) outpath=filepath
log <-  file(paste0(outpath,"_remove_duptrip_indels_23to26_log.txt"), open = "wt")
sink(log, type = "message")

# input GEN file
GEN <- fread(paste0(filepath,".gen"))
message("Input SNPs = ", nrow(GEN))

# give column names to non-GENotype columns
names(GEN)[1:5] = c("Chr","Name","MapInfo","A1","A2")
# unite Chr and MapInfo
GEN %<>% unite(Chr_MapInfo, c(Chr, MapInfo), sep="_")

# get frequency of Chr_MapInfo
gen_DupTrip <- GEN %>% count(Chr_MapInfo)
count_DupTrip = gen_DupTrip %>% count(n)
message("Unique SNPs = ", count_DupTrip[[1,2]], "*1 = ", count_DupTrip[[1,2]]*1)
message("Duplicate SNPs = ", count_DupTrip[[2,2]], "*2 = ", count_DupTrip[[2,2]]*2)
message("Triplicate SNPs = ", count_DupTrip[[3,2]], "*3 = ", count_DupTrip[[3,2]]*3)
message("Multiple SNPs = ", count_DupTrip[[2,2]]*2, "+", count_DupTrip[[3,2]]*3, " = ", count_DupTrip[[2,2]]*2+count_DupTrip[[3,2]]*3)

gen_DupTrip %>%
  ggplot(aes(x=n))+geom_bar()+geom_text(stat='count',aes(label=..count..),vjust=-0.5)+
  labs(title="Frequency of unique, duplicate, triplicate SNPs", x="Occurrence of same Chr_MapInfo", y="Number of SNPs")+theme_bw()
  ggsave(paste0(outpath,"_ChrMapInfo_frequency.png"))
  message("Output : Plot - Frequency of unique, duplicate, triplicate SNPs : ", outpath,"_ChrMapInfo_frequency.png")
# filter duplicate/triplicate Chr_MapInfo SNPs
gen_DupTrip %<>% filter(n>1) %>% select(-n) %>% inner_join(GEN, by="Chr_MapInfo")

# unique instances of duplicate/triplicate SNPs to keep
  ## select out Name because that is always unique
gen_DupTrip_keep <- gen_DupTrip %>% select(-Name) %>%
  ## only keep unique lines (Chr, MapInfo, A1, A2, GENotype) - this will retain 1 instance of identical SNPs and multiple instances of SNPs with different A1/A2/GENotype
  distinct(.keep_all=TRUE) %>%
  ## get frequency of multiple Chr_MapInfo and filter unique instances - these are the SNPs to keep
  count(Chr_MapInfo) %>% filter(n==1) %>% select(-n) %>%
  ## join with original GEN file to append names to the retained SNPs based on their Chr_MapInfo - this will add back all instances of retained SNPs because the names are different
  inner_join(GEN[,1:2], by="Chr_MapInfo") %>%
  ## only keep 1 instance of retained SNPs based on their Chr_MapInfo
  distinct(Chr_MapInfo, .keep_all=TRUE)
message("Unique instances of duplicate/triplicate (keep-list) SNPs = ", nrow(gen_DupTrip_keep))

# instances of duplicate/triplicate SNPs to discard # check below
gen_DupTrip_discard <- gen_DupTrip %>% anti_join(gen_DupTrip_keep, by="Name")
message("Redundant or non-identical instances of duplicate/triplicate (discard-list) SNPs = ", nrow(gen_DupTrip_discard))

message("Is the number of multiple SNPs equal to the sum of keep and discard SNP list? : ", count_DupTrip[[2,2]]*2+count_DupTrip[[3,2]]*3==nrow(gen_DupTrip_keep)+nrow(gen_DupTrip_discard))

# filter out instances of duplicate/triplicate SNPs to discard
gen_NODupTrip <- GEN %>% anti_join(gen_DupTrip_discard, by="Name")
message("Unique SNPs = ", nrow(gen_NODupTrip))

message("Is the number of unique SNPs equal to the difference between input and discard SNP list? : ", nrow(gen_NODupTrip)==nrow(GEN)-nrow(gen_DupTrip_discard))

# arrange final NODupTrip GEN file by Chr and then by MapInfo
gen_NODupTrip_final = gen_NODupTrip %>% arrange(Chr_MapInfo) %>% separate(Chr_MapInfo, c("Chr","MapInfo"), sep="_")
# change order of columns to reflect original order of GEN file
gen_NODupTrip_final = gen_NODupTrip_final[,c(1,3,2,4:ncol(gen_NODupTrip_final))]

write_csv(data.frame(gen_DupTrip_discard[["Name"]]), paste0(outpath, "_discard_duptrip_Names.csv"), col_names = FALSE)
message("Output : List - Names of duplicate, triplicate SNPs to discard : ", outpath,"_discard_duptrip_Names.csv")
write_csv(gen_DupTrip_discard, paste0(outpath, "_discard_duptrip_Full.csv"))
message("Output : List - Details of duplicate, triplicate SNPs to discard : ", outpath,"_discard_duptrip_Full.csv")
fwrite(gen_NODupTrip, paste0(outpath, "_NOduptrip.gen"), col.names = FALSE, row.names = FALSE, sep=" ")
message("Output : GEN - Only unique Chr and MapInfo SNPs remain : ", outpath,"_NOduptrip.gen")

# filter out SNPs with indels
gen_NODupTrip_NOindel <- gen_NODupTrip_final %>% filter(A1!="I") %>% filter(A1!="D") %>% filter(A2!="I") %>% filter(A2!="D")
# SNPs with indels only
gen_NODupTrip_indel <- gen_NODupTrip_final %>% anti_join(gen_NODupTrip_NOindel, by="Name")

message("In-del SNPs = ", nrow(gen_NODupTrip_indel))
message("Unique alleles in In-del SNPs : A1 = ", unique(unique(gen_NODupTrip_indel[,4:5])[["A1"]]), "; A2 = ", unique(unique(gen_NODupTrip_indel[,4:5])[["A2"]]))
message("Non-In-del SNPs = ", nrow(gen_NODupTrip_NOindel))
message("Unique alleles in In-del SNPs : A1 = ", unique(unique(gen_NODupTrip_NOindel[,4:5])[["A1"]]), "; A2 = ", unique(unique(gen_NODupTrip_NOindel[,4:5])[["A2"]]))

write_csv(data.frame(gen_NODupTrip_indel[["Name"]]), paste0(outpath, "_NOduptrip_discard_indel_Names.csv"), col_names = FALSE)
message("Output : List - Names of indel SNPs to discard : ", outpath,"_NOduptrip_discard_indel_Names.csv")
write_csv(gen_NODupTrip_indel, paste0(outpath, "_NOduptrip_discard_indel_Full.csv"))
message("Output : List - Details of indel SNPs to discard : ", outpath,"_NOduptrip_discard_indel_Full.csv")
fwrite(gen_NODupTrip_NOindel, paste0(outpath, "_NOduptrip_NOindel.gen"), col.names = FALSE, row.names = FALSE, sep=" ")
message("Output : GEN - Only non-indel SNPs remain : ", outpath,"_NOduptrip_NOindel.gen")

# filter out Chr>22
gen_NODupTrip_NOindel_Chr1to22 <- gen_NODupTrip_NOindel %>% mutate(Chr = as.numeric(Chr)) %>%  filter(Chr<23)
# Chr>22 only
gen_NODupTrip_NOindel_Chr23to26 <- gen_NODupTrip_NOindel %>% anti_join(gen_NODupTrip_NOindel_Chr1to22, by="Name")

message("Chr 23 to 26 SNPs = ", nrow(gen_NODupTrip_NOindel_Chr23to26))
message("Chr 1 to 22 SNPs = ", nrow(gen_NODupTrip_NOindel_Chr1to22))

write_csv(data.frame(gen_NODupTrip_NOindel_Chr23to26[["Name"]]), paste0(outpath, "_NOduptrip_NOindel_discard_Chr23to26_Names.csv"), col_names = FALSE)
message("Output : List - Names of Chr 23 to 26 SNPs to discard : ", outpath,"_NOduptrip_NOindel_discard_Chr23to26_Names.csv")
write_csv(gen_NODupTrip_NOindel_Chr23to26, paste0(outpath, "_NOduptrip_NOindel_discard_Chr23to26_Full.csv"))
message("Output : List - Details of Chr 23 to 26 SNPs to discard : ", outpath,"_NOduptrip_NOindel_discard_Chr23to26_Full.csv")
fwrite(gen_NODupTrip_NOindel_Chr1to22, paste0(outpath, "_NOduptrip_NOindel_Chr1to22.gen"), col.names = FALSE, row.names = FALSE, sep=" ")
message("Output : GEN - Only Chr 1 to 22 SNPs remain : ", outpath,"_NOduptrip_NOindel_Chr1to22.gen")

closeAllConnections()
}

remove_duptrip_indels_23to26(args[1],args[2])