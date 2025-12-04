library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)
library(ieugwasr)
library(MRPRESSO)

# Set paths and binaries
plink_bin <- genetics.binaRies::get_plink_binary()
bfile_ref <- "EUR"  # Reference genotype panel for LD clumping

# Load outcome files (lipids)
outcome_files <- list.files("AD_MR", pattern = "\\.txt$", full.names = TRUE)

# Load AD exposure GWAS
AD <- fread("AD_2022_Grch37.txt")

# Validate required columns for AD GWAS
required_cols_AD <- c("variant_id", "beta", "standard_error", "effect_allele", "other_allele",
                      "effect_allele_frequency", "p_value", "n_cases")

if (!all(required_cols_AD %in% colnames(AD))) {
  stop("AD GWAS missing required columns")
}

# Format AD as exposure
AD_exposure <- format_data(
  as.data.frame(AD),
  type = "exposure",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value",
  samplesize_col = "n_cases",
  phenotype_col = "Alzheimer's Disease"
)

# LD clumping
clumped <- tryCatch({
  ld_clump(
    dplyr::tibble(rsid = AD_exposure$SNP, pval = AD_exposure$pval.exposure),
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
  stop("No SNPs left after clumping AD exposure. Aborting.")
}

# Filter AD exposure to clumped SNPs
AD_exposure <- AD_exposure[AD_exposure$SNP %in% clumped$rsid, ]

# Calculate F-statistics for instruments
AD_exposure$F_statistic <- (AD_exposure$beta.exposure^2) / (AD_exposure$se.exposure^2)

# Initialize result lists
mr_results_list <- list()
mr_results_excl_list <- list()
pleio_results_list <- list()
hetero_results_list <- list()
egger_intercepts_list <- list()
mr_presso_results_list <- list()
mr_presso_details_list <- list()

# Loop through each lipid outcome file
for (file in outcome_files) {
  cat("\nProcessing outcome file:", file, "\n")
  
  lipid_raw <- fread(file)
  lipid_name <- gsub(".*BIO_MTB_(.*?)\\.chr.*", "\\1", basename(file))
  lipid_raw$Lipid <- lipid_name

  # Check for required outcome columns
  required_cols <- c("SNP", "ES", "SE", "A1", "A0", "MAF", "P", "N")
  if (!all(required_cols %in% colnames(lipid_raw))) {
    cat("Missing required columns in outcome file. Skipping.\n")
    next
  }

  # Match SNPs with AD exposure
  outcome_sub <- lipid_raw[SNP %in% AD_exposure$SNP, ]
  if (nrow(outcome_sub) == 0) {
    cat("No overlapping SNPs in outcome for", lipid_name, ". Skipping.\n")
    next
  }

  # Format outcome data
  outcome_dat <- format_data(
    as.data.frame(outcome_sub),
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

  # Harmonize exposure and outcome
  harmonized_dat <- tryCatch(
    harmonise_data(AD_exposure, outcome_dat),
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
  harmonized_dat$exposure <- "AD"

  # Perform MR
  mr_res <- tryCatch({
    mr(harmonized_dat)
  }, error = function(e) {
    cat("MR failed:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(mr_res)) next

  mr_res$trait <- lipid_name
  mr_res$nsnp <- nrow(harmonized_dat)
  mr_results_list[[lipid_name]] <- mr_res

  # MR-Egger intercept
  egger_res <- tryCatch(mr_egger_regression(harmonized_dat), error = function(e) NULL)
  egger_intercepts_list[[lipid_name]] <- data.frame(
    trait = lipid_name,
    intercept_egger = ifelse(is.null(egger_res), NA, egger_res$intercept),
    se_egger = ifelse(is.null(egger_res), NA, egger_res$intercept_se),
    pval_egger = ifelse(is.null(egger_res), NA, egger_res$intercept_pval)
  )

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

  # MR-PRESSO (requires ≥3 SNPs)
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
      mr_presso_results_list[[lipid_name]] <- presso_res$`Main MR results`

      if (!is.null(presso_res$`MR results`$`Outlier corrected`)) {
        mr_results_excl_list[[lipid_name]] <- presso_res$`MR results`$`Outlier corrected` %>%
          mutate(trait = lipid_name)
      }

      if (!is.null(presso_res$`MR-PRESSO results`$`Global Test`) &&
          !is.null(presso_res$`MR-PRESSO results`$`Global Test`$Pvalue) &&
          !is.null(presso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue)) {

        outlier_count <- if (!is.null(presso_res$`MR-PRESSO results`$`Outliers`)) {
          length(presso_res$`MR-PRESSO results`$`Outliers`)
        } else {
          0
        }

        mr_presso_details_list[[lipid_name]] <- data.frame(
          trait = lipid_name,
          presso_global_pval = presso_res$`MR-PRESSO results`$`Global Test`$Pvalue,
          outlier_count = outlier_count,
          distortion_pval = presso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue
        )
      } else {
        cat("MR-PRESSO: Missing test results for", lipid_name, "- Skipping detail summary.\n")
      }
    }
  } else {
    cat("Not enough SNPs for MR-PRESSO (need ≥3). Skipping.\n")
  }

  cat("Completed:", lipid_name, "\n")
}

# Combine results
combined_results <- if (length(mr_results_list) > 0) do.call(rbind, mr_results_list) else NULL
combined_results_excl <- if (length(mr_results_excl_list) > 0) do.call(rbind, mr_results_excl_list) else NULL
combined_pleio <- if (length(pleio_results_list) > 0) do.call(rbind, pleio_results_list) else NULL
combined_hetero <- if (length(hetero_results_list) > 0) do.call(rbind, hetero_results_list) else NULL
combined_egger <- if (length(egger_intercepts_list) > 0) do.call(rbind, egger_intercepts_list) else NULL
combined_presso_details <- if (length(mr_presso_details_list) > 0) do.call(rbind, mr_presso_details_list) else NULL

# Format pleiotropy results
pleio_df <- combined_pleio %>%
  select(trait, method, egger_intercept, se, pval) %>%
  rename(
    pleio_egger_intercept = egger_intercept,
    pleio_se = se,
    pleio_pval = pval
  )

# Merge all main results
merged <- combined_results %>%
  left_join(pleio_df, by = c("trait", "method")) %>%
  left_join(combined_hetero %>% select(trait, method, Q, Q_df, Q_pval), by = c("trait", "method"))

# Save outputs
fwrite(merged, "final_combined_reverse_mr_results_AD_to_lipids.csv", sep = ",")

if (!is.null(combined_results_excl)) {
  fwrite(combined_results_excl, "final_combined_reverse_mr_results_AD_to_lipids_excl_outliers.csv", sep = ",")
}

if (!is.null(combined_presso_details)) {
  fwrite(combined_presso_details, "final_combined_mr_presso_details_AD_to_lipids.csv", sep = ",")
}

cat("All done!\n")
