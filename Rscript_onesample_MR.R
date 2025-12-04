#!/usr/bin/env Rscript
options(bitmapType='cairo')
args = commandArgs(trailingOnly=TRUE)
INPUT <- args[1]
dat1<-read.table("universe_SYMBOL_ID_one_sample_MR_data_lipid_genes_snp_corrected.txt",header=T)
library(OneSampleMR)
require(glue)
library(dplyr)
o <-INPUT%>% glue
fit1 <- ivreg::ivreg(o, data = dat1)
results_dat<-as.data.frame(t(c(summary(fit1)$coefficients[2],summary(fit1)$coefficients[14],summary(fit1)$coefficients[38],summary(fit1)$diagnostics[1,][3],summary(fit1)$diagnostics[1,][4],summary(fit1)$diagnostics[2,][4],summary(fit1)$diagnostics[3,][4],summary(fit1)$r.squared)))
names(results_dat)<-c("beta","se","pvalue","F_stat","weak_instument_pval","Wu-Hausman_pval","Sargan_pval","R_squared")
results_dat$INPUT<-INPUT
write.table(results_dat,file=paste0(gsub(INPUT,pattern="[~]|[+]|[|]",replacement="_"),".corrected.txt"),row.names=F,col.names=T,sep='\t',quote=F)
