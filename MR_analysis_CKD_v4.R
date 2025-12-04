library(ieugwasr)
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)
library(RadialMR)
library(MRPRESSO)

plink_bin <- genetics.binaRies::get_plink_binary()
bfile_ref <- "EUR"  # Reference panel for LD clumping

CKD <- fread("CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt", header = TRUE)
exposure_files <- list.files("CKD_MR", pattern = "\\.txt$", full.names = TRUE)

mr_results_list <- list()
pleio_results_list <- list()
hetero_results_list <- list()
presso_results_list <- list()
fstat_summary_list <- list()

for (file in exposure_files) {
  cat("\nProcessing:", file, "\n")
  
  lipid_raw <- fread(file)
  lipid_name <- gsub(".*BIO_MTB_(.*?)\\.chr.*", "\\1", basename(file))
  lipid_raw$Lipid <- lipid_name
  
  required_cols <- c("SNP", "ES", "SE", "A1", "A0", "MAF", "P", "N")
  if (!all(required_cols %in% colnames(lipid_raw))) {
    cat("Missing required columns. Skipping.\n")
    next
  }
  
  lipid_gw <- lipid_raw[P < 5e-8]
  if (nrow(lipid_gw) == 0) {
    cat("No genome-wide significant SNPs. Skipping.\n")
    next
  }

  # --- Compute F-statistics ---
  lipid_gw$F_stat <- (lipid_gw$ES)^2 / (lipid_gw$SE)^2
  mean_F <- mean(lipid_gw$F_stat, na.rm = TRUE)
  fstat_summary_list[[lipid_name]] <- data.frame(trait = lipid_name, mean_F_stat = mean_F)
  
  exposure_dat <- format_data(
    as.data.frame(lipid_gw),
    type = "exposure",
    snp_col = "SNP",
    beta_col = "ES",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A0",
    eaf_col = "MAF",
    pval_col = "P",
    samplesize_col = "N",
    phenotype_col = "Lipid"
  )
  
  clumped <- tryCatch({
    ld_clump(
      dplyr::tibble(rsid = exposure_dat$SNP, pval = exposure_dat$pval.exposure),
      plink_bin = plink_bin,
      bfile = bfile_ref,
      clump_kb = 10000,
      clump_r2 = 0.01,
      clump_p = 5e-8
    )
  }, error = function(e) {
    cat("Clumping failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(clumped) || nrow(clumped) == 0) {
    cat("No SNPs left after clumping. Skipping.\n")
    next
  }
  
  exposure_dat <- exposure_dat[exposure_dat$SNP %in% clumped$rsid, ]
  CKD_sub <- CKD[RSID %in% exposure_dat$SNP, ]
  if (nrow(CKD_sub) == 0) {
    cat("No overlapping SNPs with CKD data. Skipping.\n")
    next
  }
  
  outcome_dat <- format_data(
    as.data.frame(CKD_sub),
    type = "outcome",
    snp_col = "RSID",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "P-value",
    samplesize_col = "n_total_sum",
    phenotype_col = "Chronic Kidney Disease"
  )
  
  harmonized <- tryCatch(harmonise_data(exposure_dat, outcome_dat), error = function(e) NULL)
  if (is.null(harmonized) || nrow(harmonized) == 0) {
    cat("Harmonisation failed or empty. Skipping.\n")
    next
  }
  
  harmonized$outcome <- "CKD"
  
  mr_res <- tryCatch(mr(harmonized), error = function(e) NULL)
  if (is.null(mr_res) || nrow(mr_res) == 0) {
    cat("MR analysis failed or returned no results. Skipping.\n")
    next
  }
  mr_res$trait <- lipid_name
  mr_res$nsnp <- nrow(harmonized)
  mr_results_list[[lipid_name]] <- mr_res
  
  pleio <- tryCatch(mr_pleiotropy_test(harmonized), error = function(e) NULL)
  if (!is.null(pleio) && nrow(pleio) > 0) {
    pleio$trait <- lipid_name
    pleio$method <- "MR Egger"
    pleio_results_list[[lipid_name]] <- pleio
  } else {
    cat("No pleiotropy results for", lipid_name, "\n")
  }
  
  hetero <- tryCatch(mr_heterogeneity(harmonized), error = function(e) NULL)
  if (!is.null(hetero) && nrow(hetero) > 0) {
    hetero$trait <- lipid_name
    hetero_results_list[[lipid_name]] <- hetero
  } else {
    cat("No heterogeneity results for", lipid_name, "\n")
  }
  
  if (nrow(harmonized) >= 4) {
    cat("Running MR-PRESSO for:", lipid_name, "\n")
    tryCatch({
      presso_res <- mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = harmonized,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      )
      
      presso_results_list[[lipid_name]] <- data.frame(
        trait = lipid_name,
        method = "MR-PRESSO",
        presso_global_pval = tryCatch(presso_res$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA),
        outlier_count = if (!is.null(presso_res$`MR-PRESSO results`$`Outliers`)) nrow(presso_res$`MR-PRESSO results`$`Outliers`) else 0,
        distortion_pval = tryCatch(presso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue, error = function(e) NA)
      )
    }, error = function(e) {
      cat("MR-PRESSO failed:", conditionMessage(e), "\n")
    })
  } else {
    cat("Not enough SNPs for MR-PRESSO in:", lipid_name, "\n")
  }
  
  cat("Completed:", lipid_name, "\n")
}

# Combine results
combined_mr <- if (length(mr_results_list) > 0) do.call(rbind, mr_results_list) else NULL
combined_pleio <- if (length(pleio_results_list) > 0) do.call(rbind, pleio_results_list) else NULL
combined_hetero <- if (length(hetero_results_list) > 0) do.call(rbind, hetero_results_list) else NULL
combined_presso <- if (length(presso_results_list) > 0) do.call(rbind, presso_results_list) else NULL
combined_fstat <- if (length(fstat_summary_list) > 0) do.call(rbind, fstat_summary_list) else NULL

# Clean and merge results
if (!is.null(combined_pleio)) {
  pleio_df <- combined_pleio %>%
    select(trait, method, egger_intercept, se, pval) %>%
    rename(
      pleio_egger_intercept = egger_intercept,
      pleio_se = se,
      pleio_pval = pval
    )
} else pleio_df <- NULL

if (!is.null(combined_hetero)) {
  hetero_df <- combined_hetero %>%
    select(trait, method, Q, Q_df, Q_pval)
} else hetero_df <- NULL

final_results <- combined_mr
if (!is.null(pleio_df)) final_results <- left_join(final_results, pleio_df, by = c("trait", "method"))
if (!is.null(hetero_df)) final_results <- left_join(final_results, hetero_df, by = c("trait", "method"))
if (!is.null(combined_presso)) final_results <- left_join(final_results, combined_presso, by = c("trait", "method"))
if (!is.null(combined_fstat)) final_results <- left_join(final_results, combined_fstat, by = "trait")

write.csv(final_results, "final_combined_mr_results_CKD_with_PRESSO.csv", row.names = FALSE)
cat("✅ All done!\n")
