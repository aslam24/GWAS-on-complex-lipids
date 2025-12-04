library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)
library(ieugwasr)
library(MRPRESSO)

plink_bin <- genetics.binaRies::get_plink_binary()
bfile_ref <- "EUR"  # Reference panel for LD clumping

exposure_files <- list.files("CAD_MR", pattern = "\\.txt$", full.names = TRUE)
CAD <- fread("METAGWAS2.tbl.rsid.chrbp")  # CAD GWAS file

mr_results_list <- list()
pleio_results_list <- list()
hetero_results_list <- list()
egger_intercepts_list <- list()
mrpresso_results_list <- list()
fstat_summary_list <- list()  # <<-- For mean F-statistics

for (file in exposure_files) {
  cat("\nProcessing file:", file, "\n")
  
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
    cat("LD clumping failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(clumped) || nrow(clumped) == 0) {
    cat("No SNPs left after clumping. Skipping.\n")
    next
  }
  
  exposure_dat <- exposure_dat[exposure_dat$SNP %in% clumped$rsid, ]
  
  CAD_sub <- CAD[oldID %in% exposure_dat$SNP, ]
  if (nrow(CAD_sub) == 0) {
    cat("No overlapping SNPs in CAD outcome. Skipping.\n")
    next
  }
  
  outcome_dat <- format_data(
    as.data.frame(CAD_sub),
    type = "outcome",
    snp_col = "oldID",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "P-value",
    samplesize_col = "N",
    phenotype_col = "CAD"
  )
  
  harmonized_dat <- tryCatch(
    harmonise_data(exposure_dat, outcome_dat),
    error = function(e) {
      cat("Harmonisation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(harmonized_dat) || !is.data.frame(harmonized_dat) || nrow(harmonized_dat) == 0) {
    cat("Invalid or empty harmonized dataset. Skipping.\n")
    next
  }
  
  harmonized_dat$outcome <- "CAD"
  nsnp <- nrow(harmonized_dat)
  
  method_list <- if (nsnp == 1) {
    "mr_wald_ratio"
  } else if (nsnp == 2) {
    c("mr_ivw", "mr_wald_ratio")
  } else {
    c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
  }
  
  mr_res <- tryCatch({
    mr(harmonized_dat, method_list = method_list)
  }, error = function(e) {
    cat("MR failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(mr_res)) next
  
  mr_res$trait <- lipid_name
  mr_res$nsnp <- nsnp
  mr_results_list[[lipid_name]] <- mr_res
  
  egger_res <- tryCatch({
    if (nsnp >= 3) mr_egger_regression(harmonized_dat) else NULL
  }, error = function(e) NULL)
  
  egger_intercepts_list[[lipid_name]] <- data.frame(
    trait = lipid_name,
    method = "MR Egger",
    intercept_egger = ifelse(is.null(egger_res), NA, egger_res$intercept),
    se_egger = ifelse(is.null(egger_res), NA, egger_res$intercept_se),
    pval_egger = ifelse(is.null(egger_res), NA, egger_res$intercept_pval)
  )
  
  pleio_res <- tryCatch({
    if (nsnp >= 3) {
      pleio <- mr_pleiotropy_test(harmonized_dat)
      pleio$trait <- lipid_name
      pleio$method <- "MR Egger"
      pleio
    } else NULL
  }, error = function(e) NULL)
  if (!is.null(pleio_res)) pleio_results_list[[lipid_name]] <- pleio_res
  
  hetero_res <- tryCatch({
    if (nsnp >= 2) {
      hetero <- mr_heterogeneity(harmonized_dat)
      hetero$trait <- lipid_name
      hetero
    } else NULL
  }, error = function(e) NULL)
  if (!is.null(hetero_res)) hetero_results_list[[lipid_name]] <- hetero_res
  
  ## --- MR-PRESSO ---
  if (nsnp >= 4) {
    cat("Running MR-PRESSO for:", lipid_name, "\n")
    tryCatch({
      presso_res <- mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = harmonized_dat,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      )
      
      global_pval <- tryCatch(presso_res$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA)
      distortion_pval <- tryCatch(presso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue, error = function(e) NA)
      outlier_count <- if (!is.null(presso_res$`MR-PRESSO results`$`Outliers`)) {
        nrow(presso_res$`MR-PRESSO results`$`Outliers`)
      } else {
        0
      }
      
      mrpresso_results_list[[lipid_name]] <- data.frame(
        trait = lipid_name,
        method = "MR-PRESSO",
        presso_global_pval = global_pval,
        outlier_count = outlier_count,
        distortion_pval = distortion_pval
      )
    }, error = function(e) {
      cat("MR-PRESSO failed:", conditionMessage(e), "\n")
    })
  }
  
  cat("Completed:", lipid_name, "\n")
}

# Combine results
combined_results <- do.call(rbind, mr_results_list)
combined_pleio <- do.call(rbind, pleio_results_list)
combined_hetero <- do.call(rbind, hetero_results_list)
combined_egger <- do.call(rbind, egger_intercepts_list)
combined_presso <- do.call(rbind, mrpresso_results_list)
combined_fstat <- do.call(rbind, fstat_summary_list)  # <<-- F-statistics

# Merge into final results
final_combined <- combined_results %>%
  left_join(
    combined_pleio %>% select(trait, method, egger_intercept, se, pval) %>%
      rename(pleio_egger_intercept = egger_intercept, pleio_se = se, pleio_pval = pval),
    by = c("trait", "method")
  ) %>%
  left_join(
    combined_hetero %>% select(trait, method, Q, Q_df, Q_pval),
    by = c("trait", "method")
  ) %>%
  left_join(
    combined_presso,
    by = c("trait", "method")
  ) %>%
  left_join(
    combined_fstat,
    by = "trait"
  )

# Export final results
fwrite(final_combined, "final_combined_mr_results_CAD_with_presso.csv", row.names = FALSE)
cat("✅ All done. Results saved.\n")
