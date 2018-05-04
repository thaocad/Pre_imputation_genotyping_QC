#Script Name: STEP2 QC pipeline
#Description: Automate all QC steps on output files of Genome Studio to create input files for Marie's pipeline.
#							This script must be accompanied with 5 extra R scripts: STEP2_plots.R, STEP8_VisualInspect.R, STEP12_Sexcheck_and_IdentitybyDescentcheck.R, STEP17_matchTOP_KIDS.R, STEP19_remove_duptrip_indels_23to26.R
#Output:1) statistical plots (for step 2, step 8 and step 12) for users to determine cut-offs
#				2) files for Marie's pipeline
#				3) Log file for QC pipeline
#Version: 1.0
#Author: Thi Thu Thao Nguyen

function error_exit
{
#	----------------------------------------------------------------
#	Function for exit due to fatal program error
#		Accepts 1 argument:
#			string containing descriptive error message
#	----------------------------------------------------------------
	echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
	exit 1
}

DATE=`date +%d%b%Y` # date of QC
### Theoretically, you need to change only the 2 following variables for the whole QC pipeline
INPUT="EXAMPLE_08Mar2018" # name of ped and map file (from Genome_Studio plink output, without file extension; for example, EXAMPLE_03Mar2018 for the 2 input files EXAMPLE_03Mar2018.map and EXAMPLE_03Mar2018.ped)
STEM="EXAMPLE_09Mar2018" # expected name of your output (without file extension)

###################################################################################################################
echo "--------STEP1-Raw_output_to_bfile--------"
mkdir -p STEP1-Raw_output_to_bfile
plink --file INPUT/$INPUT \
--make-bed \
--noweb \
--out STEP1-Raw_output_to_bfile/$STEM

###log
echo "--------STEP1-Raw_output_to_bfile--------" >>QC_$DATE.log
grep "Total genotyping rate" ./STEP1-Raw_output_to_bfile/*.log >>QC_$DATE.log
grep "variants.*people pass filters and QC." ./STEP1-Raw_output_to_bfile/*.log >> QC_$DATE.log


###################################################################################################################
echo "--------STEP2-Generate_sample_and_SNP_summary_stats--------"
mkdir -p STEP2-Generate_sample_and_SNP_summary_stats
plink --bfile STEP1-Raw_output_to_bfile/$STEM \
--missing \
--het \
--noweb \
--out STEP2-Generate_sample_and_SNP_summary_stats/$STEM.stats

sample_missing=$STEM".stats.imiss"
SNP_missing=$STEM".stats.lmiss"
het=$STEM".stats.het"

mkdir -p STEP2-Generate_sample_and_SNP_summary_stats/PLOTS
cd STEP2-Generate_sample_and_SNP_summary_stats

## running R script to generate plots to decide missing SNP/Sample rate and heterozygosity rate cut-off values
Rscript ../STEP2_plots.R $sample_missing $SNP_missing $het

#####INPUT FROM THE PLOTS
read -p "Based on SampleMissing.png plot in STEP2 PLOTS folder, what is your sample missing data rate cut-off? " mind
read -p "Based on SNPMissing.png plot in STEP2 PLOTS folder, what is your SNP missing data rate cut-off? " geno


mindtoprint=`echo $mind*100/1 | bc ` ###divided by 1 to get rid of trailing 0s
echo "sample missing-rate cut-off used in plink command is " $mind
mind1=`echo 1-$mind | bc`
echo "sample missing-rate cut-off is "$mind1


genotoprint=`echo $geno*100/1 | bc ` ###divided by 1 to get rid of trailing 0s
echo "SNP missing-rate cut-off used in plink command is " $geno
geno1=`echo 1-$geno | bc`
echo "SNP missing-rate cut-off is "$geno1

cd ..

###log
echo "" >> QC_$DATE.log
echo "--------STEP2-Generate_sample_and_SNP_summary_stats--------" >>QC_$DATE.log
echo "Chosen Sample Missing Data Rate Cut-Off was "$mind1 >> QC_$DATE.log
echo "Chosen Variant-based(SNP) Missing Data Rate Cut-Off was "$geno1 >> QC_$DATE.log

####################################################################################################################
echo "--------STEP3-Remove_samples_with_greater_than_"$mindtoprint"_percent_missing_data--------"
mkdir -p "STEP3-Remove_samples_with_greater_than_"$mindtoprint"_percent_missing_data"
plink --bfile STEP1-Raw_output_to_bfile/$STEM \
--mind $mind \
--make-bed \
--noweb \
--out "STEP3-Remove_samples_with_greater_than_"$mindtoprint"_percent_missing_data/"$STEM".GSAMPLE"

###log
echo "" >> QC_$DATE.log
echo "--------STEP3-Remove_samples_with_greater_than_"$mindtoprint"_percent_missing_data--------" >>QC_$DATE.log
grep "Total genotyping rate" STEP3-Remove_samples_with_greater_than_"$mindtoprint"_percent_missing_data/*.log >>QC_$DATE.log
grep "people removed due to missing genotype data" STEP3-Remove_samples_with_greater_than_"$mindtoprint"_percent_missing_data/*.log >> QC_$DATE.log
grep "variants.*people pass filters and QC." STEP3-Remove_samples_with_greater_than_"$mindtoprint"_percent_missing_data/*.log >> QC_$DATE.log

####################################################################################################################
echo "--------STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data--------"
mkdir -p "STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data"
plink --bfile "STEP3-Remove_samples_with_greater_than_"$mindtoprint"_percent_missing_data/"$STEM".GSAMPLE" \
--geno $geno \
--make-bed \
--noweb \
--out "STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data/"$STEM".GSNP"

###log
echo "" >> QC_$DATE.log
echo "--------STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data--------" >>QC_$DATE.log
grep "Total genotyping rate" STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data/*.log >>QC_$DATE.log
grep "variants removed due to missing genotype data" STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data/*.log >> QC_$DATE.log
grep "variants.*people pass filters and QC." STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data/*.log >> QC_$DATE.log


#######################################################################################################################
echo "--------STEP 5-6-7 of MAVAN: SKIP HERE, NO NEED TO MERGE--------"
mkdir -p STEP_5-6-7_EMPTY
printf "\n"

###log
echo "" >> QC_$DATE.log
echo "--------STEP 5-6-7 of MAVAN: SKIP HERE, NO NEED TO MERGE--------" >>QC_$DATE.log

#########################################################################################################################
echo "--------STEP8-Analyze_and_Visualize_Summary_Statistics--------"
mkdir -p STEP8-Analyze_and_Visualize_Summary_Statistics
plink \
--bfile STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data/"$STEM".GSNP \
--missing \
--freq \
--hardy \
--het \
--noweb \
--out STEP8-Analyze_and_Visualize_Summary_Statistics/$STEM.stats

mkdir -p STEP8-Analyze_and_Visualize_Summary_Statistics/PLOTS
cd STEP8-Analyze_and_Visualize_Summary_Statistics


hwe=$STEM".stats.hwe"
frq=$STEM".stats.frq"

Rscript ../STEP8_VisualInspect.R $sample_missing $SNP_missing $hwe $het $frq

cd ..

if [ "$?" != "0" ]; then
  error_exit "***** Error ******"
fi

# Some decisions were made based on visual inspection of the stats output
read -p "Based on SNP_MAF.png plot in STEP8 PLOTS folder, what is your MAF cut-off? " maf_cutoff
read -p "Based on SNP_HWE.png plot in STEP2 PLOTS folder, what is your HWE cut-off? " hwe_cutoff
hwe_cutoff_toprint=$hwe_cutoff
hwe_cutoff="1e-"$hwe_cutoff
echo $hwe_cutoff_toprint
echo $hwe_cutoff

###log
echo "" >> QC_$DATE.log
echo "--------STEP8-Analyze_and_Visualize_Summary_Statistics--------" >>QC_$DATE.log
echo "Chosen MAF cut-off was "$maf_cutoff >> QC_$DATE.log
echo "Chosen HWE cut-off was "$hwe_cutoff >> QC_$DATE.log


####################################################################################################################
# STEP9
# Purpose: To remove SNPs with low MAF and odd HWE
printf "\n"
echo "--------STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE--------"
mkdir -p STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE
plink \
--bfile STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data/"$STEM".GSNP \
--maf $maf_cutoff \
--hwe $hwe_cutoff \
--make-bed \
--noweb \
--out STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE/$STEM"_maf_hwe"

###log
echo "" >> QC_$DATE.log
echo "--------STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE--------" >>QC_$DATE.log
grep "Total genotyping rate" STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE/*.log  >> QC_$DATE.log
grep "variants removed due to Hardy-Weinberg exact test" STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE/*.log >> QC_$DATE.log
grep "variants removed due to minor allele" STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE/*.log >> QC_$DATE.log
grep "variants and.*pass filters and QC" STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE/*.log >> QC_$DATE.log

if [ "$?" != "0" ]; then
  error_exit "***** Error ******"
fi


########################################################################################################################
printf "\n"
echo "--------STEP10-LD_pruning--------"
# to generate SNP list to be kept
mkdir -p STEP10-LD_pruning
plink \
--bfile STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE/$STEM"_maf_hwe" \
--indep-pairwise 50 5 0.5 \
--noweb \
--out STEP10-LD_pruning/$STEM"_maf_hwe.indeppair"

###log
echo "" >> QC_$DATE.log
echo "--------STEP10-LD_pruning--------" >>QC_$DATE.log
grep "Total genotyping rate" STEP10-LD_pruning/"$STEM"_maf_hwe.indeppair.log  >> QC_$DATE.log
grep "Pruned.*variants from chromosome" STEP10-LD_pruning/"$STEM"_maf_hwe.indeppair.log >> QC_$DATE.log
grep "variants removed" STEP10-LD_pruning/"$STEM"_maf_hwe.indeppair.log >> QC_$DATE.log

if [ "$?" != "0" ]; then
  error_exit "***** Error ******"
fi

# LD-pruning
plink \
--bfile STEP9-Remove_SNPs_with_low_MAF_and_odd_HWE/$STEM"_maf_hwe" \
--extract STEP10-LD_pruning/$STEM"_maf_hwe.indeppair.prune.in" \
--make-bed \
--out STEP10-LD_pruning/$STEM"_maf_hwe_prune"
if [ "$?" != "0" ]; then
  error_exit "***** Error ******"
fi

###log_continue
grep "Total genotyping rate" STEP10-LD_pruning/"$STEM"_maf_hwe_prune.log >> QC_$DATE.log
grep "variants and.*pass filters and QC" STEP10-LD_pruning/"$STEM"_maf_hwe_prune.log >> QC_$DATE.log


########################################################################################################################
# STEP11
# Purpose: To produce gender and identity by descent results from available genotyping data

printf "\n"
echo "--------STEP11-Sex_check_and_Identity_by_descent_check--------"
## SEX-CHECK
mkdir -p STEP11-Sex_check_and_Identity_by_descent_check
mkdir -p STEP11-Sex_check_and_Identity_by_descent_check/_SexCheck
mkdir -p STEP11-Sex_check_and_Identity_by_descent_check/_IdentitybyDescent

plink \
--bfile STEP10-LD_pruning/$STEM"_maf_hwe_prune" \
--check-sex \
--noweb \
--out STEP11-Sex_check_and_Identity_by_descent_check/_SexCheck/$STEM"_maf_hwe_prune.sex"
if [ "$?" != "0" ]; then
  error_exit "***** Error ******"
fi
## IDENTITY-BY-DESCENT CHECK
plink \
--bfile STEP10-LD_pruning/$STEM"_maf_hwe_prune" \
--genome \
--noweb \
--out STEP11-Sex_check_and_Identity_by_descent_check/_IdentitybyDescent/$STEM"_maf_hwe_prune.genome"

if [ "$?" != "0" ]; then
  error_exit "***** Error ******"
fi

###log
echo "" >> QC_$DATE.log
echo "--------STEP11-Sex_check_and_Identity_by_descent_check--------" >>QC_$DATE.log


########################################################################################################################
# STEP12,
#Purpose: This is a sanity check for samples to confirm they are identified correctly until this step.
#Method: Mislabeled samples were identified and reported in this step.

printf "\n"
echo "--------STEP12-Generate_a_clean_list_of_samples--------"
mkdir -p STEP12-Generate_a_clean_list_of_samples
read -p "What is the name of the file containing gender information (.txt, stored in INPUT folder, with 2 columns ID and Gender)? " gender
read -p "What is the value for female in gender file (example: 0,1, or 2)? " female
read -p "What is the value for male in gender file (example: 0,1, or 2)? " male
read -p "What is the name of the file containing sample information (.csv, original sample sheet to link Sentrix ID and Sample ID, stored in INPUT folder)? " sample

imiss=STEP8-Analyze_and_Visualize_Summary_Statistics/"$STEM".stats.imiss
genome=STEP11-Sex_check_and_Identity_by_descent_check/_IdentitybyDescent/"$STEM"_maf_hwe_prune.genome.genome
sexcheck=STEP11-Sex_check_and_Identity_by_descent_check/_SexCheck/"$STEM"_maf_hwe_prune.sex.sexcheck
Rscript STEP12_Sexcheck_and_IdentitybyDescentcheck.R $imiss $genome $sexcheck $gender $sample $female $male


########################################################################################################################
# STEP13
# Purpose: Some decisions should be made based on the misspecified SNPs (sex/IdentitybyDescent)
# Method: Mislabeled samples (if exist) idendified from STEP12 will be excluded in this step
printf "\n"
printf "\n"
printf "\n"
read -p "Based on data from the txt files on Sex and Identity by Descent in STEP12 folder, \
please make a text file in STEP12 folder listing of all samples which should be removed, \
including their_UniqueID [tab] reason_to_be_removed(optional) , one \"UniqueID\"\t\"reason\" per line, then name the file as SampleExclude_clean.txt. \
For example:
9479477167_R10C02 misspecified_female_and_duplicate
9479476111_R04C01 misspecified_female
9479621052_R05C02 unmatched_male

Note for SexCheck:
  F score for FEMALE is less that 0.2 and annotated value (CHILDSEX) should be 1, and F score for MALE is more than 0.8 and annotated value (CHILDSEX) should be 0
  \"unmatch\" means CHILDSEX and F score are not reflecting the same gender or CHILDSEX is missing, \"misspecified\" means unmatch and oppositely annotated (female was annotated as male and vice versa)
  \"uncertain\" means F score is in an uncertain zone (0.2<=F<=0.8)

Note for Identity_by_Descent
  Duplicate_or_identical_twin.txt will notify users about the duplicate samples or potential identical twins, among which only one sample should be retained.

Do you have any troubled samples, such as misspecified sex or surprised duplicate samples, to be removed \
Is the SampleExclude_clean.txt file ready? answer NO when there is no sample to be excluded?" yes_or_no

while [[ ! "$yes_or_no" =~ ^([yY][eE][sS]|[yY]|[nN][oO]|[nN])+$ ]]
do
read -p "Wrong answer! Is the SampleExclude_clean.txt file ready? answer NO when there is no sample to be excluded?" yes_or_no
done

case $yes_or_no in
###case 1: there are troubled samles to be removed
  [yY][eE][sS]|[yY])
    tr -cd '\11\12\40-\176' < STEP12-Generate_a_clean_list_of_samples/SampleExclude_clean.txt > temp
    awk 'NR==FNR{a[$1];next} $2 in a {print $1,$2}' temp STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data/"$STEM".GSNP.fam > STEP12-Generate_a_clean_list_of_samples/SampleExclude.txt

    #log
    echo "" >> QC_$DATE.log
    echo "--------STEP12-Generate_a_clean_list_of_samples--------" >>QC_$DATE.log
    echo "The following samples were excluded:">>  QC_$DATE.log
    cat temp >>  QC_$DATE.log
    echo "" >> QC_$DATE.log
    rm temp*



    printf "\n"
    echo "--------STEP13-Remove_questionable_samples--------"
    mkdir -p STEP13-Remove_questionable_samples
    plink --bfile STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data/"$STEM".GSNP \
    --remove STEP12-Generate_a_clean_list_of_samples/SampleExclude.txt \
    --noweb \
    --make-bed \
    --out STEP13-Remove_questionable_samples/$STEM"_samprem"

    output_step13=STEP13-Remove_questionable_samples/$STEM"_samprem"

    if [ "$?" != "0" ]; then
      error_exit "***** Error ******"
    fi

    #log
    echo "" >> QC_$DATE.log
    echo "--------STEP13-Remove_questionable_samples--------" >>QC_$DATE.log
    grep "-remove:" STEP13-Remove_questionable_samples/$STEM"_samprem.log"  >> QC_$DATE.log
    grep "Total genotyping rate" STEP13-Remove_questionable_samples/$STEM"_samprem.log"  >> QC_$DATE.log
    grep "variants and.*pass filters and QC" STEP13-Remove_questionable_samples/$STEM"_samprem.log"  >> QC_$DATE.log

    ##recheck statistics in STEP13
    cd STEP13-Remove_questionable_samples
    mkdir -p PLOTS
    plink --bfile "$STEM"_samprem --missing --freq --hardy --het --noweb --out "$STEM"_samprem.stats
    sample_missing="$STEM"_samprem.stats.imiss
    SNP_missing="$STEM"_samprem.stats.lmiss
    hwe="$STEM"_samprem.stats.hwe
    het="$STEM"_samprem.stats.het
    frq="$STEM"_samprem.stats.frq
    Rscript ../STEP8_VisualInspect.r $sample_missing $SNP_missing $hwe $het $frq


    read -p "Please check PLOTS in step13, then enter the new hwe cut off?" hwe_cutoff
    hwe_cutoff_toprint=$hwe_cutoff
    hwe_cutoff="1e-"$hwe_cutoff

    cd ..
    #log
    echo "" >> QC_$DATE.log
    echo "Chosen HWE cut-off after removing troubled samples was "$hwe_cutoff >> QC_$DATE.log
  ;;


###case 2: there are no troubled samples to be removed, STEP13 will be skipped
  [nN][oO]|[nN])
    echo "--------STEP13:SKIPPED, no troubled samples--------"
    mkdir -p STEP_13_EMPTY
    printf "\n"
    echo "" >> QC_$DATE.log
    echo "--------STEP13:SKIPPED, no troubled samples--------" >>QC_$DATE.log

    output_step13=STEP4-Remove_SNPS_with_greater_than_"$genotoprint"_percent_missing_data/"$STEM".GSNP
  ;;
esac

########################################################################################################################
printf "\n"
echo "--------STEP14-HWE"$hwe_cutoff_toprint"--------"
mkdir -p STEP14-HWE"$hwe_cutoff_toprint"
plink --bfile $output_step13 \
--hwe $hwe_cutoff \
--make-bed \
--noweb \
--out "STEP14-HWE"$hwe_cutoff_toprint"/"$STEM"_samprem_"$hwe_cutoff_toprint
if [ "$?" != "0" ]; then
  error_exit "***** Error ******"
fi

#log
echo "" >> QC_$DATE.log
echo "--------STEP14-HWE"$hwe_cutoff_toprint"--------" >>QC_$DATE.log
grep "hwe:" STEP14-HWE"$hwe_cutoff_toprint"/"$STEM"_samprem_"$hwe_cutoff_toprint".log  >> QC_$DATE.log
grep "variants and.*pass filters and QC" STEP14-HWE"$hwe_cutoff_toprint"/"$STEM"_samprem_"$hwe_cutoff_toprint".log  >> QC_$DATE.log


########################################################################################################################
read -p "Was any chr 0 in your array   \
(If the answer is yes, then please have the file snp0bp_removelist.txt ready in the INPUT folder.)?" yes_or_no_1

while [[ ! "$yes_or_no_1" =~ ^([yY][eE][sS]|[yY]|[nN][oO]|[nN])+$ ]]
do
read -p "Wrong answer! Was PsychChip array (or any other array rather than PsychArray) used for this cohort?" yes_or_no_1
done

case $yes_or_no_1 in

###case 1: PsychChip array was used, so no need to remove chr 0
  [nN][oO]|[nN])

  printf "\n"
  mkdir -p STEP_15_EMPTY
  echo "-----------STEP15-Remove_chr_0: SKIPPED, reason: no chromosome 0 SNPs in the SNP list--------------"
  output_step15="STEP14-HWE"$hwe_cutoff_toprint"/"$STEM"_samprem_"$hwe_cutoff_toprint

  #log
  echo "" >> QC_$DATE.log
  echo "--------STEP15-Remove_chr_0: SKIPPED, reason: no chromosome 0 SNPs in the SNP list--------" >>QC_$DATE.log
  ;;


###case 2: PsychArray has been used
  [yY][eE][sS]|[yY])

  printf "\n"
  echo " --------STEP15-Remove_chr_0------"
  mkdir -p STEP15-Remove_chr_0
  ###remove SNPs with unknown chromosome (chr 0)
  plink --bfile "STEP14-HWE"$hwe_cutoff_toprint"/"$STEM"_samprem_"$hwe_cutoff_toprint \
  --exclude INPUT/snp0bp_removelist.txt \
  --noweb \
  --make-bed \
  --out STEP15-Remove_chr_0/$STEM"_samprem_hwe40-0bp"

  output_step15=STEP15-Remove_chr_0/$STEM"_samprem_hwe40-0bp"

    if [ "$?" != "0" ]; then
    error_exit "*****Error ******"
    fi
  #log

  ;;
esac



########################################################################################################################
printf "\n"

read -p "Based on data from the txt files on potential siblings in STEP12 folder, please make a text file in STEP12 folder listing of all samples which should be removed, including their UniqueID (Sentrix ID), one UniqueID per line, named Remove_Siblings.txt. \

Is there any potential sibling in this study?" yes_or_no_2


while [[ ! "$yes_or_no_2" =~ ^([yY][eE][sS]|[yY]|[nN][oO]|[nN])+$ ]]
do
read -p "Wrong answer! Is there any potential sibling in this study?" yes_or_no_1
done

case $yes_or_no_2 in
  ###case 1: there are potential siblings to be removed
  [yY][eE][sS]|[yY])
  echo " --------STEP16_NOSIB------"
  mkdir -p STEP16_NOSIB

  tr -cd '\11\12\40-\176' < STEP12-Generate_a_clean_list_of_samples/Remove_Siblings.txt > temp
  paste -d ' ' temp temp > STEP12-Generate_a_clean_list_of_samples/Remove_Siblings.txt
  rm temp


  plink --bfile $output_step15 \
  --exclude STEP12-Generate_a_clean_list_of_samples/Remove_Siblings.txt \
  --noweb \
  --make-bed \
  --out STEP16_NOSIB/$STEM"_samprem_hwe"$hwe_cutoff_toprint"-0bp_nosib"

  output_step16=STEP16_NOSIB/$STEM"_samprem_hwe"$hwe_cutoff_toprint"-0bp_nosib"

  if [ "$?" != "0" ]; then
    error_exit "*****Error ******"
  fi
  echo "------------------SHOULD ADD SOME DOCUMENTATION BASED ON LOG FILES------------------------"
  ;;

  ###case 2: there are no siblings to be removed
  [nN][oO]|[nN])
  mkdir -p STEP16_SKIPPED_NOSIB_to_be_removed
  printf "\n"
  echo "------------------ STEP16_SKIPPED_NOSIB_to_be_removed ------------------------"
  output_step16=$output_step15

  #log
  echo "" >> QC_$DATE.log
  echo "--------STEP16_SKIPPED_NOSIB_to_be_removed--------" >>QC_$DATE.log
  ;;
esac
########################################################################################################################
printf "\n"
echo "--------STEP17-TOP--------"
mkdir -p STEP17-TOP



plink --bfile $output_step16 \
--recode oxford \
--out STEP17-TOP/"$STEM"_step17


plink --bfile $output_step16 \
--recode \
--out STEP17-TOP/"$STEM"_step17

mapfile=STEP17-TOP/"$STEM"_step17.map
genfile=STEP17-TOP/"$STEM"_step17.gen
read -p "What is the name of the file containing update_alleles information (.update_alleles.txt, stored in RESOURCES folder, usually obtained from Will Rayner website)? " topstrand
topstrandfile=RESOURCES/$topstrand

Rscript STEP17_matchTOP_KIDS.R $mapfile $genfile $topstrandfile

plink --data STEP17-TOP/"$STEM"_step17 --a1-allele STEP17-TOP/Array_TOP1.txt --recode oxford --make-bed --out STEP17-TOP/"$STEM"_step17_TOP1
plink --data STEP17-TOP/"$STEM"_step17_TOP1 --a1-allele STEP17-TOP/Array_TOP2.txt --recode oxford --make-bed --out STEP17-TOP/"$STEM"_step17_TOP1_TOP2


echo "" >> QC_$DATE.log
echo "--------STEP17-TOP--------" >>QC_$DATE.log

########################################################################################################################
printf "\n"
echo "--------STEP18_NOMULTI_SKIPPED--------"
mkdir -p STEP18_NOMULTI_SKIPPED
output_step18=STEP17-TOP/"$STEM"_step17_TOP1_TOP2

echo "" >> QC_$DATE.log
echo "--------STEP18_NOMULTI_SKIPPED--------" >>QC_$DATE.log

########################################################################################################################
printf "\n"
echo "--------STEP19_RMDUP--------"
mkdir -p STEP19_RMDUP

output_step19=STEP19_RMDUP/"$STEM"_step19
Rscript STEP19_remove_duptrip_indels_23to26.R $output_step18 $output_step19

echo "" >> QC_$DATE.log
echo "--------STEP19_RMDUP--------" >>QC_$DATE.log
cat STEP19_RMDUP/*_log.txt >> QC_$DATE.log

plink \
--bfile $output_step18 \
--exclude STEP19_RMDUP/"$STEM"_step19_discard_duptrip_Names.csv \
--make-bed \
--keep-allele-order \
--out STEP19_RMDUP/"$STEM"_step19_nomulti_nodup


plink \
--bfile STEP19_RMDUP/"$STEM"_step19_nomulti_nodup \
--exclude STEP19_RMDUP/"$STEM"_step19_NOduptrip_discard_indel_Names.csv \
--make-bed \
--keep-allele-order \
--out STEP19_RMDUP/"$STEM"_step19_nomulti_nodup_noindel

plink \
--bfile STEP19_RMDUP/"$STEM"_step19_nomulti_nodup_noindel \
--exclude STEP19_RMDUP/"$STEM"_step19_NOduptrip_NOindel_discard_Chr23to26_Names.csv \
--recode \
--keep-allele-order \
--out STEP19_RMDUP/"$STEM"_imputePrepInput_"$DATE"


echo "" >> QC_$DATE.log
echo "--------QC HAS BEEN SUCCESSFULLY FINISHED--------" >>QC_$DATE.log

echo "--------QC HAS BEEN SUCCESSFULLY FINISHED--------"
