library(ieugwasr)
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)
library(RadialMR)
library(MRPRESSO)

plink_bin <- genetics.binaRies::get_plink_binary()
bfile_ref <- "EUR"  # Reference panel for LD clumping

# Load outcome data
T2D <- fread("T2D_Xue_et_al_2018.txt")

# List exposure GWAS summary stats files
exposure_files <- list.files("T2D_MR", pattern = "\\.txt$", full.names = TRUE)

# Containers for results
mr_results_list <- list()
mr_results_excl_list <- list()
pleio_results_list <- list()
pleio_results_excl_list <- list()
hetero_results_list <- list()
hetero_results_excl_list <- list()
egger_intercepts_list <- list()
mr_presso_results_list <- list()

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
  
  T2D_sub <- T2D[SNP %in% exposure_dat$SNP, ]
  if (nrow(T2D_sub) == 0) {
    cat("No overlapping SNPs in T2D outcome. Skipping.\n")
    next
  }
  
  outcome_dat <- format_data(
    as.data.frame(T2D_sub),
    type = "outcome",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "frq_A1",
    pval_col = "P",
    samplesize_col = "N",
    phenotype_col = "Type 2 Diabetes"
  )
  
  harmonized_dat <- tryCatch(
    harmonise_data(exposure_dat, outcome_dat),
    error = function(e) {
      cat("Harmonisation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(harmonized_dat) || nrow(harmonized_dat) == 0) {
    cat("Invalid or empty harmonized dataset. Skipping.\n")
    next
  }
  
  harmonized_dat$outcome <- "T2D"
  harmonized_dat$F_statistic <- (harmonized_dat$beta.exposure^2) / (harmonized_dat$se.exposure^2)
  mean_F <- mean(harmonized_dat$F_statistic, na.rm = TRUE)
  
  mr_res <- tryCatch({
    mr(harmonized_dat)
  }, error = function(e) {
    cat("MR failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(mr_res) || nrow(mr_res) == 0) {
    cat("No MR results. Skipping.\n")
    next
  }
  
  mr_res$trait <- lipid_name
  mr_res$nsnp <- nrow(harmonized_dat)
  mr_res$F_statistic <- mean_F
  mr_results_list[[lipid_name]] <- mr_res
  
  pleio_res <- tryCatch(mr_pleiotropy_test(harmonized_dat), error = function(e) NULL)
  if (!is.null(pleio_res) && nrow(pleio_res) > 0) {
    pleio_res$trait <- lipid_name
    pleio_res$method <- "MR Egger"
    pleio_results_list[[lipid_name]] <- pleio_res
  }
  
  hetero_res <- tryCatch(mr_heterogeneity(harmonized_dat), error = function(e) NULL)
  if (!is.null(hetero_res) && nrow(hetero_res) > 0) {
    hetero_res$trait <- lipid_name
    hetero_results_list[[lipid_name]] <- hetero_res
  }
  
  egger_res <- tryCatch(mr_egger_regression(harmonized_dat), error = function(e) NULL)
  if (!is.null(egger_res) && nrow(egger_res) > 0) {
    egger_intercepts_list[[lipid_name]] <- data.frame(
      trait = lipid_name,
      intercept_egger = egger_res$intercept,
      se_egger = egger_res$intercept_se,
      pval_egger = egger_res$intercept_pval
    )
  }
  
  # MR-PRESSO
  if (nrow(harmonized_dat) >= 3) {
    presso_res <- tryCatch({
      mr_presso(BetaOutcome = "beta.outcome",
                BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome",
                SdExposure = "se.exposure",
                OUTLIERtest = TRUE,
                DISTORTIONtest = TRUE,
                data = harmonized_dat,
                NbDistribution = 1000,
                SignifThreshold = 0.05)
    }, error = function(e) {
      cat("MR-PRESSO failed:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (!is.null(presso_res)) {
      global_p <- ifelse(!is.null(presso_res$`MR-PRESSO results`$`Global Test`$Pvalue),
                         presso_res$`MR-PRESSO results`$`Global Test`$Pvalue, NA)
      distortion_p <- ifelse(!is.null(presso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue),
                             presso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue, NA)
      outlier_n <- ifelse(!is.null(presso_res$`MR-PRESSO results`$`Outliers`),
                          nrow(presso_res$`MR-PRESSO results`$`Outliers`), 0)
      
      mr_presso_results_list[[lipid_name]] <- data.frame(
        trait = lipid_name,
        presso_global_pval = global_p,
        distortion_pval = distortion_p,
        outlier_count = outlier_n
      )
      
      if (!is.null(presso_res$`MR results`$`Outlier corrected`)) {
        excl_mr <- presso_res$`MR results`$`Outlier corrected` %>%
          mutate(trait = lipid_name, F_statistic = mean_F)
        mr_results_excl_list[[lipid_name]] <- excl_mr
      }
    }
  }
  
  cat("Completed:", lipid_name, "\n")
}

# Combine all
combined_results <- if (length(mr_results_list) > 0) do.call(rbind, mr_results_list) else NULL
combined_results_excl <- if (length(mr_results_excl_list) > 0) do.call(rbind, mr_results_excl_list) else NULL
combined_pleio <- if (length(pleio_results_list) > 0) do.call(rbind, pleio_results_list) else NULL
combined_pleio_excl <- if (length(pleio_results_excl_list) > 0) do.call(rbind, pleio_results_excl_list) else NULL
combined_hetero <- if (length(hetero_results_list) > 0) do.call(rbind, hetero_results_list) else NULL
combined_hetero_excl <- if (length(hetero_results_excl_list) > 0) do.call(rbind, hetero_results_excl_list) else NULL
combined_presso <- if (length(mr_presso_results_list) > 0) do.call(rbind, mr_presso_results_list) else NULL

combined_pleio_all <- rbind(na.omit(combined_pleio), na.omit(combined_pleio_excl))
combined_hetero_all <- rbind(na.omit(combined_hetero), na.omit(combined_hetero_excl))

pleio_df <- combined_pleio_all %>%
  select(trait, method, egger_intercept, se, pval) %>%
  rename(pleio_egger_intercept = egger_intercept, pleio_se = se, pleio_pval = pval)

hetero_df <- combined_hetero_all %>%
  select(trait, method, Q, Q_df, Q_pval)

if (!is.null(combined_results)) {
  final_combined <- combined_results %>%
    left_join(pleio_df, by = c("trait", "method")) %>%
    left_join(hetero_df, by = c("trait", "method")) %>%
    left_join(combined_presso, by = "trait")
  
  write.csv(final_combined, "final_combined_mr_results_T2D.csv", row.names = FALSE)
}

if (!is.null(combined_results_excl)) {
  final_combined_excl <- combined_results_excl %>%
    left_join(pleio_df, by = c("trait", "method")) %>%
    left_join(hetero_df, by = c("trait", "method")) %>%
    left_join(combined_presso, by = "trait")
  
  write.csv(final_combined_excl, "final_combined_mr_results_excl_outliers_T2D.csv", row.names = FALSE)
}

cat("All done!\n")
