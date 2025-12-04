library(ieugwasr)
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)
library(RadialMR)
library(MRPRESSO)

plink_bin <- genetics.binaRies::get_plink_binary()
bfile_ref <- "EUR"  # Reference panel for clumping

# Load CKD data as exposure
CKD <- fread("CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt", header = TRUE)
CKD_gw <- CKD[`P-value` < 5e-8]

# Format CKD exposure data
exposure_dat <- format_data(
  as.data.frame(CKD_gw),
  type = "exposure",
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P-value",
  samplesize_col = "n_total_sum",
  phenotype_col = "CKD"
)

# LD clumping CKD SNPs
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
  cat("LD clumping failed for CKD:", conditionMessage(e), "\n")
  return(NULL)
})

if (is.null(clumped) || nrow(clumped) == 0) {
  stop("No SNPs left after clumping CKD exposure")
}
exposure_dat <- exposure_dat[exposure_dat$SNP %in% clumped$rsid, ]

# Load lipid outcome files
outcome_files <- list.files("CAD_MR", pattern = "\\.txt$", full.names = TRUE)

# Initialize result lists
mr_results_list <- list()
pleio_results_list <- list()
hetero_results_list <- list()
egger_intercepts_list <- list()
presso_results_list <- list()

for (file in outcome_files) {
  cat("\nProcessing outcome file:", file, "\n")
  
  lipid_raw <- fread(file)
  lipid_name <- gsub(".*BIO_MTB_(.*?)\\.chr.*", "\\1", basename(file))
  
  required_cols <- c("SNP", "ES", "SE", "A1", "A0", "MAF", "P", "N")
  if (!all(required_cols %in% colnames(lipid_raw))) {
    cat("Missing required columns. Skipping.\n")
    next
  }
  
  lipid_outcome <- lipid_raw[SNP %in% exposure_dat$SNP]
  if (nrow(lipid_outcome) == 0) {
    cat("No overlapping SNPs in outcome for", lipid_name, "\n")
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
  
  harmonized_dat <- tryCatch(
    harmonise_data(exposure_dat, outcome_dat),
    error = function(e) {
      cat("Harmonisation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(harmonized_dat) || nrow(harmonized_dat) == 0) {
    cat("No harmonised data. Skipping.\n")
    next
  }
  
  harmonized_dat$exposure <- "CKD"
  harmonized_dat$outcome <- lipid_name
  
  mr_res <- tryCatch(mr(harmonized_dat), error = function(e) NULL)
  if (!is.null(mr_res) && nrow(mr_res) > 0) {
    mr_res$trait <- lipid_name
    mr_res$nsnp <- nrow(harmonized_dat)
    mr_results_list[[lipid_name]] <- mr_res
  } else {
    cat("MR failed for", lipid_name, "\n")
    next
  }
  
  # Egger intercept
  egger_res <- tryCatch(mr_egger_regression(harmonized_dat), error = function(e) NULL)
  if (!is.null(egger_res)) {
    egger_intercepts_list[[lipid_name]] <- data.frame(
      trait = lipid_name,
      intercept_egger = egger_res$intercept,
      se_egger = egger_res$intercept_se,
      pval_egger = egger_res$intercept_pval
    )
  }
  
  # Pleiotropy test
  pleio_res <- tryCatch(mr_pleiotropy_test(harmonized_dat), error = function(e) NULL)
  if (!is.null(pleio_res) && nrow(pleio_res) > 0) {
    pleio_res$trait <- lipid_name
    pleio_res$method <- "MR Egger"
    pleio_results_list[[lipid_name]] <- pleio_res
  }
  
  # Heterogeneity test
  hetero_res <- tryCatch(mr_heterogeneity(harmonized_dat), error = function(e) NULL)
  if (!is.null(hetero_res) && nrow(hetero_res) > 0) {
    hetero_res$trait <- lipid_name
    hetero_results_list[[lipid_name]] <- hetero_res
  }

  # MR-PRESSO
  if (nrow(harmonized_dat) >= 4) {
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
      if (!is.null(presso_res$`Main MR results`)) {
        presso_results_list[[lipid_name]] <- data.frame(
          trait = lipid_name,
          presso_causal_estimate = presso_res$`Main MR results`$`Causal Estimate`,
          presso_raw_pval = presso_res$`Main MR results`$`P-value`,
          presso_distortion_pval = presso_res$`Distortion Test`$`P-value Distortion`
        )
      }
    }, error = function(e) {
      cat("MR-PRESSO failed for", lipid_name, ":", conditionMessage(e), "\n")
    })
  }

  cat("Completed:", lipid_name, "\n")
}

# Combine and export
combined_results <- if(length(mr_results_list) > 0) do.call(rbind, mr_results_list) else NULL
combined_pleio <- bind_rows(pleio_results_list)
combined_hetero <- bind_rows(hetero_results_list)
combined_egger <- bind_rows(egger_intercepts_list)
combined_presso <- bind_rows(presso_results_list)

pleio_df <- if (!is.null(combined_pleio)) {
  combined_pleio %>%
    select(trait, method, egger_intercept, se, pval) %>%
    rename(
      pleio_egger_intercept = egger_intercept,
      pleio_se = se,
      pleio_pval = pval
    )
} else NULL

hetero_df <- if (!is.null(combined_hetero)) {
  combined_hetero %>%
    select(trait, method, Q, Q_df, Q_pval)
} else NULL

final_combined <- combined_results
if (!is.null(combined_egger)) final_combined <- left_join(final_combined, combined_egger, by = "trait")
if (!is.null(pleio_df)) final_combined <- left_join(final_combined, pleio_df, by = c("trait", "method"))
if (!is.null(hetero_df)) final_combined <- left_join(final_combined, hetero_df, by = c("trait", "method"))
if (!is.null(combined_presso)) final_combined <- left_join(final_combined, combined_presso, by = "trait")

if (!is.null(final_combined)) {
  write.csv(final_combined, "reverse_mr_results_CKD_to_lipids_with_presso.csv", row.names = FALSE)
}

cat("Reverse MR with MR-PRESSO complete!\n")
