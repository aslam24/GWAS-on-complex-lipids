#!/bin/bash
#####i## Pipeline to Run GWAS on complex lipids
MAX_JOBS=200

launch()
{
    while [ $(jobs | wc -l) -ge $MAX_JOBS ]
    do
        # at least $MAX_JOBS are still running.
        sleep 2 
    done
    "$@" &
}

WORKDIR="/groups/rs-projects/GCL"
# Diretory to indicate the Log files from GWAS to be stored.
logdirectory="${WORKDIR}/log_III_1_12"
# Directory to merge chromosome 1-22 for each metabolite and save it to this folder.
resultsdirectory="${WORKDIR}/results_III_1_12"
# Directory to store metabolite-phenotype names list , snp list and Actual phenotype list
datadirectory="${WORKDIR}/data_2"
# Directory to store Genotype Vcf file.
dosedirectory="${WORKDIR}/dose_2"

study="RS1"
ethnic="EA"
platform="Nightingale"
date="20250226"
name="Imtiaz"
phenotypefile="RS_EUR_Grch37_Model3_nobatch.rvtest.txt"
infofile="snp.info"

cd $logdirectory
#we set the outer loop to chromosomes and the inner loop to phenotypes
for chr in {1..22};

do

  for phenotype in `awk '{print$1}' OFS='\t' ${datadirectory}/${study}.${ethnic}.${platform}.linkfile.${date}.${name}.txt`; 
   do
   ### prepare the format that we need
   echo "CHR POS EFFECT_Allele NON_Effect_ALLELE EAF BETA SE PVALUE Rsq HWE CALL_RATE"  > ${resultsdirectory}/${study}.${ethnic}.${platform}.${phenotype}.${date}.${name}.txt

      echo "`date`: rvtest:	$phenotype		$chr." >> RunGWAS_rs1.log
      launch rvtest \
       --inVcf ${dosedirectory}/Omni_final_ver2_updated_${chr}.phased.impute2_final.vcf.gz \
       --pheno ${datadirectory}/${phenotypefile} \
       --pheno-name ${phenotype} \
       --dosage DS \
       --siteFile ${datadirectory}/Omni_final_ver2_updated_${chr}.phased.impute2_final.snplist \
       --siteMACMin 5 \
       --out ${logdirectory}/${study}.${ethnic}.${platform}.${phenotype}.chr${chr}.${date}.${name} \
       --meta score
   done
   wait
done
####### Merging chromosome 1 to 22 for each phenotype ############
#for phenotype in `awk '{print$1}' OFS='\t' ${datadirectory}/${study}.${ethnic}.${platform}.linkfile.${date}.${name}.txt`;
#do
#  cat ${logdirectory}/${study}.${ethnic}.${platform}.${phenotype}.chr*.${date}.${name}.MetaScore.assoc | grep -v "^#" | grep -v "^CHROM" > ${logdirectory}/${study}.${ethnic}.${platform}.${phenotype}.chr1_22.${date}.${name}.MetaScore.assoc
#done
