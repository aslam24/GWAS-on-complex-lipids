#!/usr/bin/env Rscript
options(bitmapType='cairo')
args = commandArgs(trailingOnly=TRUE)
INPUT <- args[1]
exposure <- strsplit(INPUT,":")[[1]][1]
outcome <- strsplit(INPUT,":")[[1]][2]
#print(strsplit(INPUT,":")[[1]][1])
#print(strsplit(INPUT,":")[[1]][2])
print(exposure)
print(outcome)
library(TwoSampleMR)
library(ieugwasr)
library(genetics.binaRies)
library(dplyr)
#options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
library(data.table)
exposure_data<-read.table(paste0(exposure),header=T)
exposure_dat <- read_exposure_data(
  filename = paste0(exposure),
  sep = '\t',
  snp_col = 'SNP',
  beta_col = 'slope',
  se_col = 'slope_se',
  effect_allele_col = 'A1',
  phenotype_col = 'SYMBOL',
  units_col = 'units',
  other_allele_col = 'A2',
  eaf_col = 'maf',
  samplesize_col = 'ma_samples',
  pval_col = 'pval_nominal'
)
print(head(exposure_dat))
exposure_dat<-exposure_dat[exposure_dat$pval.exposure<1e-4,]
#options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
#exposure_dat_1 <- clump_data(exposure_dat)
Clumped_data<-ld_clump(dplyr::tibble(rsid=exposure_dat$SNP,pval=exposure_dat$pval.exposure),plink_bin = genetics.binaRies::get_plink_binary(),bfile="EUR",clump_kb = 10000,clump_r2 = 0.01,clump_p= 1e-4)
exposure_dat<-exposure_dat[exposure_dat$SNP%in%Clumped_data$rsid,]
exposure_dat$Fstat<- exposure_dat$beta.exposure^2 / exposure_dat$se.exposure^2
meanF<-mean(exposure_dat$Fstat)
write.table(meanF,file=paste0(INPUT,".","F_stat_2025.txt"),row.names=F,col.names=F,sep='\t',quote=F)

#CHROM	POS	REF	ALT	N_INFORMATIVE	AF	INFORMATIVE_ALT_AC	CALL_RATE	HWE_PVALUE	N_REF	N_HET	N_ALT	U_STAT	SQRT_V_STAT	ALT_EFFSIZE	PVALUE
#outcome_data <-fread(paste0(outcome),header = T)
outcome_data <-fread(paste("grep -v '^#'", outcome),header=F)
names(outcome_data)<-c("CHROM"	,"POS"	,"REF"	,"ALT"	,"N_INFORMATIVE"	,"AF"	,"INFORMATIVE_ALT_AC"	,"CALL_RATE"	,"HWE_PVALUE"	,"N_REF"	,"N_HET"	,"N_ALT"	,"U_STAT"	,"SQRT_V_STAT"	,"ALT_EFFSIZE"	,"PVALUE")
outcome_data$outcome <- outcome
outcome_data$phenotype <- outcome
outcome_data$SNPID<-paste0(outcome_data$CHROM,":",outcome_data$POS)
outcome_data$SE<-1/outcome_data$SQRT_V_STAT
snp_data<-exposure_data[,c(1,25)]
outcome_data<-outcome_data[outcome_data$SNPID%in%snp_data$SNPID,]
outcome_data<-merge(outcome_data,snp_data,by="SNPID")
outcome_data<-as.data.frame(outcome_data)
outcome_data_gwas <- format_data(
  outcome_data,
  type="outcome",
  snp_col = "SNP",
  beta_col = "ALT_EFFSIZE",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "PVALUE",
  chr_col = "CHROM",
  pos_col = "POS",
  phenotype_col = "phenotype")
outcome_data_gwas$outcome<-gsub(outcome_data_gwas$outcome,pattern = "RS1.EA.Nightingale.",replacement = "")
outcome_data_gwas$outcome<-gsub(outcome_data_gwas$outcome,pattern = ".20250226.Imtiaz.MetaScore.assoc",replacement = "")
dat <- harmonise_data(exposure_dat, outcome_data_gwas, action = 2)

write.table(dat,file=paste0(exposure,"_",outcome,".","harmonized_lipid_species_2025.txt"),row.names=F,col.names=T,sep='\t',quote=F)
mr_results <- mr(dat)
print(mr_results)
mr_results$F<-meanF
write.table(mr_results,file=paste0(exposure,"_",outcome,".","mr_lipid_species_2025.txt"),row.names=F,col.names=T,sep='\t',quote=F)
