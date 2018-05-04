This document provides detailed instructions for STEP2-QC and STEP3-PRE_IMPUTATION (Marie’s pipeline) to prepare files ready for imputation from plink data (from GENOME_STUDIO processed data). Only files for STEP2-QC are provided in this repository.   
# STEP2-QC    
## All following scripts should be placed in the STEP2-QC folder:
-	runQC.sh: main running script
-	supported R scripts: STEP2_plots.R, STEP8_VisualInspect.R, STEP12_Sexcheck_and_IdentitybyDescentcheck.R, STEP17_matchTOP_KIDS.R, STEP19_remove_duptrip_indels_23to26.R    
## Required input: 
1)	In INPUT subfolder
-	.map and .ped files from Genome Studio
-	Gender data: text file with 2 columns, ID and gender (0/1 or 1/2) 
-	Sample sheet (.csv)
-	SNP list with chr 0 (snp0bp_removelist.txt) to be removed
2)	In RESOURCE subfolder:
-	Update alleles file from Will Rayner website (http://www.well.ox.ac.uk/~wrayner/strand/ABtoTOPstrand.html)    
## Running process:
1)	You should change the 2 variables at the beginning of the main running script (runQC.sh): INPUT and STEM (using Vim, Nano or other linux command line editor; in case you use a Window editor such as Notepad or Atom, you have to run ‘dos2unix runQC.sh’ before running the script)
2)	In Linux command line window, type: bash runQC.sh
3)	The script will be paused at several places, so you can check the plots and decide the value cut-offs for filtering parameters such as maf, missing call rate, etc.   
a.	First pause is at STEP 2: You need to check the plots in STEP2 PLOTs folder to determine the cut-off values for Sample Missing Rate and SNP Missing Rate then enter these directly in the command line window (for example, Sample Missing Rate cut-off is 0.01 for 99% Sample Call Rate and SNP Missing Rate cut-off is 0.05 for 95% SNP Call Rate    
b.	Second pause is at STEP 8: You need to check the plots in STEP8 PLOTS folder to determine the cut-off values for MAF and HWE (Normally, MAF cut-off is at 0.05, and HWE cut-off is usually determined to remove outliers in the cohort)   
c.	Third pause is at STEP 12: You need to provide the name of your gender file (stored in INPUT folder, for example, EXAMPLE_genderdata.txt), the value for female in your file (0,1 or 2), the value of male in your file (0,1 or 2), and name of your sample sheet file (for example, EXAMPLE_Samplesheet.csv)   
d.	Fourth pause is at STEP 13: You need to check in STEP 12 folder for mis-annotated samples (mis-specified sex), duplicate or identical twin samples and prepare a file named SampleExclude_clean.txt in STEP 12 then answer Yes for the question asked in the command line. If there is no bad sample to be removed, answer No.   
e.	Fifth pause is at STEP 15: If there are SNPs with chromosome 0 (maybe due to conflicts with previous manifests, Illumina manifest has SNPs with chr 0), we usually try to remove them to eliminate future complications. Please prepare a file with a list of SNP names with chr 0, and named the file as snp0bp_removelist.txt (INPUT folder)   
f.	Sixth pause is at STEP 16: If there are potential siblings (you can check in the file first_degree.txt in STEP12), prepare the list of their Unique ID (Sentrix ID), one per line, named it as Remove_Siblings.txt and saved it in STEP 12 folder, answer Yes for the question, otherwise, answer No    
g.	Seventh pause is at STEP 18: Please provide the name of the update alleles file in RESOURCE folder     


# STEP3-PRE_IMPUTATION:
Marie Forest created this pipeline to format files ready for imputation (Sanger Imputation Service). Please refer to her github page for more information: https://github.com/eauforest/imputePrepSanger  
All following scripts should be placed in the STEP3-PRE_IMPUTATION folder:   
-	imputePrep_script.sh: main running script   
-	supported scripts: HRC-1000G-check-bim_v4.2.7.pl, reportRedaction.sh, update_build.sh   
-	optional scripts: imputePrepExe.sh (to record the command and flags that you used to run the main script), README.md (detailed information on how to run the script)   
Required input:    
## input subfolder:    
*	all files \*_imputePrepInput_[date].\* from step 19 of STEP2-QC   
*	.miss, . muplitple, .strand files from Will Rayner website http://www.well.ox.ac.uk/~wrayner/strand/   
## resources subfolder:    
*	HRC.r1-1.GRCh37.wgs.mac5.sites.tab (downloaded from http://www.haplotype-reference-consortium.org/site)   
*	human_g1k_v37.fasta/ human_g1k_v37.fasta.gz (downloaded from   ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz)   
*	human_g1k_v37.fasta.fai   
*	ucsc2ensembl.txt    
## Command to run the script (should be recorded in imputePrepExe.sh for future reference):     
./imputePrep_script.sh LEARS_09Mar2018_imputePrepInput_13Mar2018 \ # step 19 files, without ext.
GSA-24v1-0_A2-b37 \ # file names of .miss, . multiple, .strand files without ext.
1.0 1.0 0.05 1e-99 \ # --mind --geno --maf –hwe (we gave max value because STEP2-QC have done this)    
input  \ # input folder name    
resources \ #resources folder name   
results # expected output folder name      
    
Output file to submit to IMPUTATION SANGER SERVER is the file .vcf.gz from the folder ‘results’ in Marie’s pipeline
