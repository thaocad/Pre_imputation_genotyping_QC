library(data.table)
library(tidyverse)
library(gridExtra)
library(magrittr)
args = commandArgs(trailingOnly=TRUE)

# ##generate gender file
# Case_Control_File_for_Kieran <- read_sav("//cortex.meaney.lab/pub/Thao_LEARS/STEP2-QC/INPUT/Case Control File for Kieran.sav")
# Gender <- Case_Control_File_for_Kieran %>% select(ID, Sex)
# colnames(Gender) <- c("SampleID","gender_0male_1female")
# write_tsv(Gender,"LEARS_genderdata.txt")

imiss <- fread(args[1],  header=TRUE)

#Identity by descent
data <- fread(args[2], header=TRUE)
genderfile <- paste0("INPUT/",args[4],sep="")
samplefile <- paste0("INPUT/",args[5],sep="")

####FROM ANNOTATED DATA
gender<- fread(genderfile, header=TRUE)
colnames(gender)<- c("Study_ID","CHILDSEX")
m<- as.numeric(args[7])
f<- as.numeric(args[6])

gender[gender$CHILDSEX==m,]$CHILDSEX<-0
gender[gender$CHILDSEX==f,]$CHILDSEX<-1
sample <-read_csv(samplefile,skip=11)
names(sample)[2]<- "X2"
names(sample)[3]<- "X3"
names(sample)[5]<- "X5"
sample$ID <- paste(sample$X2,"_",sample$X3, sep="")
sample<-sample[,c("X5","ID")]
colnames(sample)<- c("Sample_ID","UniqueID")

##special for LEARS case because Sample_ID was not number
sample$Sample_ID <- gsub("PSYC00","",sample$Sample_ID)
sample$Sample_ID <- gsub("PSYCH0","",sample$Sample_ID)
sample$Sample_ID <- gsub("PSCY00","",sample$Sample_ID)
sample$Sample_ID <- as.numeric(sample$Sample_ID)

##A first-degree relative is defined as a close blood relative which includes the individual's parents, full siblings, or children
##A second-degree relative is defined as a blood relative which includes the individual's grandparents, grandchildren, aunts, uncles, nephews, nieces or half-siblings
##A third-degree relative is defined as a blood relative which includes the individual's first-cousins, great-grandparents or great grandchildren
stranger<- data[data$PI_HAT<=0.1,]
identical<-  data[data$PI_HAT>0.6,]
onedegree<- data[(data$PI_HAT>=0.4)&(data$PI_HAT<=0.6),] ### siblings/ parents, half siblings
onedegree$FID1 <- sample[sample$UniqueID%in%onedegree$IID1,]$Sample_ID
onedegree$FID2 <- sample[sample$UniqueID%in%onedegree$IID2,]$Sample_ID
twodegree<-data[((data$PI_HAT>=0.2)&(data$PI_HAT<0.4)),] ### grandparents, uncles/aunts
twodegree$FID1 <- sample[sample$UniqueID%in%twodegree$IID1,]$Sample_ID
twodegree$FID2 <- sample[sample$UniqueID%in%twodegree$IID2,]$Sample_ID
threedegree<-data[((data$PI_HAT>0.1)&(data$PI_HAT<0.2)),]### cousins
threedegree$FID1 <- sample[sample$UniqueID%in%threedegree$IID1,]$Sample_ID
threedegree$FID2 <- sample[sample$UniqueID%in%threedegree$IID2,]$Sample_ID
write.table(onedegree,"STEP12-Generate_a_clean_list_of_samples/first_degree.txt", quote=F, col.names=T, row.names=F)
write.table(twodegree,"STEP12-Generate_a_clean_list_of_samples/second_degree.txt", quote=F, col.names=T, row.names=F)
write.table(threedegree,"STEP12-Generate_a_clean_list_of_samples/third_degree.txt", quote=F, col.names=T, row.names=F)

identical1 <- merge(identical,sample, by.x="IID1", by.y="UniqueID") 
identical2 <- merge(identical,sample, by.x="IID2", by.y="UniqueID") 
identical_pair<- data.frame(identical1$IID1, identical1$Sample_ID,identical2$IID1, identical2$Sample_ID)
colnames(identical_pair)<- c("UniqueID-1","SampleID-1","UniqueID-2","SampleID-2")
write.table(identical_pair,"STEP12-Generate_a_clean_list_of_samples/Duplicate_or_identical_twin.txt", quote=F, col.names=T, row.names=F)

#SEX
sex <- read.table(args[3], as.is=TRUE, header=TRUE)
male <- sex[sex$F>0.8,]
female <- sex[sex$F<0.2,]
unsure <- sex[(sex$F>=0.2)&(sex$F<=0.8),]

sample_gender <- merge(sample,gender,by.x="Sample_ID", by.y="Study_ID")


orig_male <- sample_gender[sample_gender$CHILDSEX==0,]
orig_male <- orig_male[orig_male$UniqueID %in% sex$IID,]
male_check <- merge(orig_male, male, by.x="UniqueID",by.y="IID")

orig_female<- sample_gender[sample_gender$CHILDSEX==1,]
orig_female <- orig_female[orig_female$UniqueID %in% sex$IID,]
female_check <- merge(orig_female, female, by.x="UniqueID",by.y="IID")

if (length(unsure[,1])==0){

  mis_anno_male <- orig_male[!orig_male$UniqueID %in% male$IID,]
  mis_anno_male$F<- female[female$IID %in% mis_anno_male$UniqueID,]$F
  mis_anno_female <- orig_female[!orig_female$UniqueID %in% female$IID,]
  mis_anno_female$F<- male[male$IID %in% mis_anno_female$UniqueID,]$F
  mis_anno <- rbind (mis_anno_male, mis_anno_female)
  write.table(mis_anno, "STEP12-Generate_a_clean_list_of_samples/mis_annotated_sample.txt", quote=F, col.names=T,row.names =F, sep="\t")
} else{
  mis_anno_male <- orig_male[!orig_male$UniqueID %in% male$IID,]
  mis_anno_male$F<- female[female$IID %in% mis_anno_male$UniqueID,]$F
  
  mis_anno_male2 <- orig_male[!orig_male$UniqueID %in% unsure$IID,]
  mis_anno_male2$F<- unsure[unsure$IID %in% mis_anno_male$UniqueID,]$F
  
  mis_anno_female <- orig_female[!orig_female$UniqueID %in% female$IID,]
  mis_anno_female$F<- male[male$IID %in% mis_anno_female$UniqueID,]$F
  
  mis_anno_female2 <- orig_female[!orig_female$UniqueID %in% unsure$IID,]
  mis_anno_female2$F<- unsure[unsure$IID %in% mis_anno_female$UniqueID,]$F
  mis_anno <- rbind (mis_anno_male, mis_anno_male2,mis_anno_female,mis_anno_female2)
  write.table(mis_anno, "STEP12-Generate_a_clean_list_of_samples/mis_annotated_sample.txt", quote=F, col.names=T,row.names =F, sep="\t")
  }

