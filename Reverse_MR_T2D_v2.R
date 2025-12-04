library(ieugwasr)
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)
library(MRPRESSO)  # MR-PRESSO added

# Define PLINK binary and LD reference
plink_bin <- genetics.binaRies::get_plink_binary()
bfile_ref <- "EUR"  # Reference panel for LD clumping

# Input files
exposure_files <- list.files("T2D_MR", pattern = "\\.txt$", full.names = TRUE)
T2D <- fread("T2D_Xue_et_al_2018.txt")  # Your T2D GWAS

# Result storage
mr_results_list <- list()
pleio_results_list <- list()
hetero_results_list <- list()
egger_intercepts_list <- list()
presso_results_list <- list()

# Main loop: reverse MR (T2D → lipids)
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
  
  # Filter T2D for genome-wide significant SNPs
  T2D_exposure <- T2D[P < 5e-8]
  if (nrow(T2D_exposure) == 0) {
    cat("No genome-wide significant SNPs in T2D exposure. Skipping.\n")
    next
  }
  
  # Format T2D exposure data
  exposure_dat <- format_data(
    as.data.frame(T2D_exposure),
    type = "exposure",
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
  
  # LD clumping
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
  
  # Format lipid outcome data
  lipid_outcome <- lipid_raw[SNP %in% exposure_dat$SNP, ]
  if (nrow(lipid_outcome) == 0) {
    cat("No overlapping SNPs in lipid outcome. Skipping.\n")
    next
  }
  
  outcome_dat <- format_data(
    as.data.frame(lipid_outcome),
    type = "outcome",
    snp_col = "SNP",
    beta_col = "ES",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A0",
    eaf_col = "MAF",
    pval_col = "P",
    samplesize_col = "N",
    phenotype_col = lipid_name
  )
  
  # Harmonisation
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
  
  harmonized_dat$outcome <- lipid_name
  harmonized_dat$exposure <- "T2D"
  
  # Mendelian Randomization
  mr_res <- tryCatch({
    mr(harmonized_dat)
  }, error = function(e) {
    cat("MR failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(mr_res)) next
  
  mr_res$trait <- paste0("T2D_to_", lipid_name)
  mr_res$nsnp <- nrow(harmonized_dat)
  mr_results_list[[paste0("T2D_to_", lipid_name)]] <- mr_res
  
  # Egger Intercept
  egger_res <- tryCatch(mr_egger_regression(harmonized_dat), error = function(e) NULL)
  if (!is.null(egger_res) && nrow(egger_res) > 0) {
    egger_intercepts_list[[paste0("T2D_to_", lipid_name)]] <- data.frame(
      trait = paste0("T2D_to_", lipid_name),
      intercept_egger = egger_res$intercept,
      se_egger = egger_res$intercept_se,
      pval_egger = egger_res$intercept_pval
    )
  } else {
    egger_intercepts_list[[paste0("T2D_to_", lipid_name)]] <- data.frame(
      trait = paste0("T2D_to_", lipid_name),
      intercept_egger = NA,
      se_egger = NA,
      pval_egger = NA
    )
  }
  
  # Pleiotropy
  pleio_res <- tryCatch(mr_pleiotropy_test(harmonized_dat), error = function(e) NULL)
  if (!is.null(pleio_res) && nrow(pleio_res) > 0) {
    pleio_res$trait <- paste0("T2D_to_", lipid_name)
    pleio_res$method <- "MR Egger"
    pleio_results_list[[paste0("T2D_to_", lipid_name)]] <- pleio_res
  }
  
  # Heterogeneity
  hetero_res <- tryCatch(mr_heterogeneity(harmonized_dat), error = function(e) NULL)
  if (!is.null(hetero_res) && nrow(hetero_res) > 0) {
    hetero_res$trait <- paste0("T2D_to_", lipid_name)
    hetero_results_list[[paste0("T2D_to_", lipid_name)]] <- hetero_res
  }
  
  # MR-PRESSO
  presso_res <- tryCatch({
    run_mr_presso <- mr_presso(
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
    
    if (!is.null(run_mr_presso$`Main MR results`)) {
      data.frame(
        trait = paste0("T2D_to_", lipid_name),
        presso_causal_estimate = run_mr_presso$`Main MR results`$`Causal Estimate`,
        presso_pval_global = run_mr_presso$`MR-PRESSO results`$`Global Test`$Pvalue,
        presso_distortion_pval = run_mr_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
      )
    } else {
      data.frame(
        trait = paste0("T2D_to_", lipid_name),
        presso_causal_estimate = NA,
        presso_pval_global = NA,
        presso_distortion_pval = NA
      )
    }
  }, error = function(e) {
    data.frame(
      trait = paste0("T2D_to_", lipid_name),
      presso_causal_estimate = NA,
      presso_pval_global = NA,
      presso_distortion_pval = NA
    )
  })
  
  presso_results_list[[paste0("T2D_to_", lipid_name)]] <- presso_res
  
  cat("Completed reverse MR for:", lipid_name, "\n")
}

# Combine all result tables, check length before do.call
combined_results <- if(length(mr_results_list) > 0) do.call(rbind, mr_results_list) else NULL
combined_pleio <- if(length(pleio_results_list) > 0) do.call(rbind, pleio_results_list) else NULL
combined_hetero <- if(length(hetero_results_list) > 0) do.call(rbind, hetero_results_list) else NULL
combined_presso <- if(length(presso_results_list) > 0) do.call(rbind, presso_results_list) else NULL

# Prepare pleiotropy and heterogeneity tables safely
pleio_df <- NULL
if (!is.null(combined_pleio)) {
  required_cols_pleio <- c("trait", "method", "egger_intercept", "se", "pval")
  if (all(required_cols_pleio %in% colnames(combined_pleio))) {
    pleio_df <- combined_pleio %>%
      select(trait, method, egger_intercept, se, pval) %>%
      rename(
        pleio_egger_intercept = egger_intercept,
        pleio_se = se,
        pleio_pval = pval
      )
  }
}

hetero_df <- NULL
if (!is.null(combined_hetero)) {
  required_cols_hetero <- c("trait", "method", "Q", "Q_df", "Q_pval")
  if (all(required_cols_hetero %in% colnames(combined_hetero))) {
    hetero_df <- combined_hetero %>%
      select(trait, method, Q, Q_df, Q_pval)
  }
}

# Defensive joins: start with combined_results and add if exists
final_combined <- combined_results

if (!is.null(final_combined)) {
  if (!is.null(pleio_df) && all(c("trait", "method") %in% colnames(final_combined))) {
    final_combined <- left_join(final_combined, pleio_df, by = c("trait", "method"))
  }
  
  if (!is.null(hetero_df) && all(c("trait", "method") %in% colnames(final_combined))) {
    final_combined <- left_join(final_combined, hetero_df, by = c("trait", "method"))
  }
  
  if (!is.null(combined_presso) && "trait" %in% colnames(final_combined) && "trait" %in% colnames(combined_presso)) {
    final_combined <- left_join(final_combined, combined_presso, by = "trait")
  }
}

# Write output if data exists
if (!is.null(final_combined)) {
  write.csv(final_combined, "final_combined_reverse_mr_results_T2D_to_lipids.csv", row.names = FALSE)
  cat("All done with reverse MR including MR-PRESSO!\n")
} else {
  cat("No results to write!\n")
}
