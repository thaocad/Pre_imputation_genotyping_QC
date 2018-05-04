############### Missing Rates for Genotype QC Step I
library(data.table)
library(tidyverse)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

###Missing Samples Data
imiss <- fread(args[1],header=TRUE) #merged imiss
#imiss <- fread("LEARS_09Mar2018.stats.imiss", header = T)

png(file="./PLOTS/SampleMissing.png")
ggplot(imiss, aes(x=FID, y=F_MISS)) + 
  geom_point(alpha=0.5, size=2) + 
  theme_bw() +
  labs(x="ID", y="Sample Missing Rate") + 
  geom_hline(yintercept = 0.01, colour="red")
dev.off()

#************************

###Missing SNP data
lmiss <- fread(args[2], header=TRUE)
#lmiss <- fread("LEARS_09Mar2018.stats.lmiss", header=TRUE)

png(file="./PLOTS/SNPMissing.png", width = 600, height = 480)
ggplot(lmiss) + geom_histogram(aes(x=F_MISS), fill = "black") + 
  theme_bw() + 
  geom_vline(xintercept=0.05, color="red")+  ## change the intercept line here if needed
  labs(x="SNP Missing Rate", y="Frequency") + 
  ggtitle("SNP Missing Rate Histogram")

dev.off()

#*************************
  
###HWE data
hwe <- fread(args[3], header=TRUE)
#hwe <- fread("LEARS_09Mar2018.stats.hwe", header=TRUE)

###HET data
het <- fread(args[4], header=TRUE)
#het <- fread("LEARS_09Mar2018.stats.het", header=TRUE)

###MAF data
frq <- fread(args[5], header=TRUE)
#frq <- fread("LEARS_09Mar2018.stats.frq", header=TRUE)

###QC

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

SNPstats<- merge(frq, lmiss, by=c("CHR","SNP"))
SNPstats <- merge(SNPstats,hwe, by = c("CHR","SNP"))

SNPstats$index <- rownames(SNPstats)
SNPstats$HWE <- SNPstats$P

###SNP QC graphs
png(file="./PLOTS/SNP_HWE.png")
par(mfrow=c(1,2))
plot(SNPstats$index, -log10(SNPstats$HWE), ylim = c(0,(max(-log10(SNPstats$HWE))+5)))
hist(SNPstats$HWE)
dev.off()

png(file="./PLOTS/SNP_MAF.png")
par(mfrow=c(1,2))
hist(SNPstats$MAF, main="MAF all SNPs")
hist(SNPstats$MAF, xlim=c(0, 0.1), main="MAF zoomed in")
abline(v=0.05, col="red")
dev.off()


