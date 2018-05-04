############### Missing Rates for Genotype QC Step I 
###Toggle all folder/file names as needed
library(data.table)
library(tidyverse)
library(gridExtra)
args = commandArgs(trailingOnly=TRUE)

###Missing Samples Data 
imiss <- fread(args[1],header=TRUE)
#imiss <- fread("LEARS_09Mar2018.stats.imiss", header = T)
head(imiss); dim(imiss)

##Save plot in PLOTS folder
png(file="./PLOTS/SampleMissing.png")
  ggplot(imiss, aes(x=FID, y=F_MISS)) + 
    geom_point(alpha=0.5, size=2) + 
    theme_bw() +
    labs(x="ID", y="Sample Missing Rate") + 
    geom_hline(yintercept = 0.01, colour="red")
dev.off()

###Missing SNP data
lmiss <- fread(args[2], header=TRUE)
#lmiss <- fread("LEARS_09Mar2018.stats.lmiss", header=TRUE)

head(lmiss); dim(lmiss) 

png(file="./PLOTS/SNPMissing.png", width = 600, height = 480)
  a<- ggplot(lmiss) + geom_histogram(aes(x=F_MISS), fill = "black") + 
      theme_bw() + 
      geom_vline(xintercept=0.05, color="red")+  ## change the intercept line here if needed
      labs(x="SNP Missing Rate", y="Frequency") + 
      ggtitle("SNP Missing Rate Histogram")
  b<- ggplot(lmiss) + geom_histogram(aes(x=F_MISS), fill="black") + 
      theme_bw() + 
      labs(x="SNP Missing Rate", y="Frequency") +  
      coord_cartesian(xlim=c(0,0.1))+
      ggtitle("SNP Missing Rate Histogram - Zoom in") + 
      geom_vline(xintercept = 0.05, colour = "red") ## change the intercept line here if needed
grid.arrange(a,b,nrow =1)
dev.off()

###HET data
het <- fread(args[3], header=TRUE)
#het <- fread("LEARS_09Mar2018.stats.het", header=TRUE)

Samplestats<-merge(imiss,het, by="IID")

#heterozygosity outliers: more than 3SD away from mean.
mean_het <- mean(Samplestats$F,na.rm=T)
std_het <- sd(Samplestats$F,na.rm=T)

#tag all samples with heterozygosity rate +/- 3 standard deviations from the mean
Samplestats<- Samplestats %>% mutate(outlier = as.factor(ifelse((F>(mean_het+3*std_het)|F<(mean_het -3*std_het)), 1, 0)))

#Missing vs Heterozygosity
png(file="./PLOTS/Sample_Heterozygosity.png")

ggplot(Samplestats,aes(x=F_MISS, y=F)) +
  geom_point(size=2, aes(colour=outlier), show.legend = F) + 
  scale_colour_manual(values=c("black", "red"))+ ### red samples can be considered outliers
  geom_text(data= subset(Samplestats,outlier==1), aes(label=IID), vjust= 0, hjust =0)+
  xlab("Sample Missing Rate") + ylab("Heterozygosity")+
  theme_bw()+
  geom_hline(yintercept= (mean_het+3*std_het), color= "red")+
  geom_hline(yintercept = (mean_het-3*std_het), color="red")

dev.off()
