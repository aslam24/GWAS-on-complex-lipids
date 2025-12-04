library(ieugwasr)
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)
library(MRPRESSO)

# Setup paths and tools
plink_bin <- genetics.binaRies::get_plink_binary()
bfile_ref <- "EUR"  # Reference panel for clumping

# Load exposure GWAS (CAD)
CAD <- fread("METAGWAS2.tbl.rsid.chrbp")

# Check required columns and set sample size if missing
required_cols_exp <- c("oldID", "Effect", "StdErr", "Allele1", "Allele2", "Freq1", "P-value", "N")
if (!"N" %in% colnames(CAD)) {
  cat("Sample size column 'N' missing. Adding default value of 100000.\n")
  CAD$N <- 100000
}
if (!all(required_cols_exp %in% colnames(CAD))) {
  stop("CAD GWAS is missing required columns.")
}

# Filter significant SNPs
CAD_gw <- CAD[`P-value` < 5e-8]
if (nrow(CAD_gw) == 0) stop("No genome-wide significant SNPs found.")

# Format as exposure data
exposure_dat <- format_data(
  as.data.frame(CAD_gw),
  type = "exposure",
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
  stop("LD clumping failed: ", conditionMessage(e))
})
exposure_dat <- exposure_dat[exposure_dat$SNP %in% clumped$rsid, ]

# Load lipid outcome files
lipid_files <- list.files("CAD_MR", pattern = "\\.txt$", full.names = TRUE)

# Prepare lists for results
mr_results_list <- list()
pleio_results_list <- list()
hetero_results_list <- list()
egger_intercepts_list <- list()
presso_results_list <- list()

# Loop through outcomes sequentially
for (file in lipid_files) {
  cat("\nProcessing lipid outcome file:", file, "\n")
  
  lipid_raw <- fread(file)
  lipid_name <- gsub(".*BIO_MTB_(.*?)\\.chr.*", "\\1", basename(file))
  lipid_raw$Lipid <- lipid_name
  
  required_cols_out <- c("SNP", "ES", "SE", "A1", "A0", "MAF", "P", "N")
  if (!all(required_cols_out %in% colnames(lipid_raw))) {
    cat("Skipping file with missing columns:", file, "\n")
    next
  }
  
  lipid_sub <- lipid_raw[SNP %in% exposure_dat$SNP]
  if (nrow(lipid_sub) == 0) next

  outcome_dat <- format_data(
    as.data.frame(lipid_sub),
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

  harmonized_dat <- tryCatch(harmonise_data(exposure_dat, outcome_dat), error = function(e) NULL)
  if (is.null(harmonized_dat) || nrow(harmonized_dat) == 0) next

  harmonized_dat$outcome <- lipid_name

  # F-statistic
  harmonized_dat$F_statistic <- (harmonized_dat$beta.exposure / harmonized_dat$se.exposure)^2
  Fstat_value <- mean(harmonized_dat$F_statistic, na.rm = TRUE)

  # MR analysis
  mr_res <- tryCatch(mr(harmonized_dat), error = function(e) NULL)
  if (is.null(mr_res)) next
  mr_res$trait <- lipid_name
  mr_res$nsnp <- nrow(harmonized_dat)
  mr_res$F_statistic <- Fstat_value
  mr_results_list[[lipid_name]] <- mr_res

  # Egger intercept
  egger_res <- tryCatch(mr_egger_regression(harmonized_dat), error = function(e) NULL)
  egger_intercepts_list[[lipid_name]] <- data.frame(
    trait = lipid_name,
    intercept_egger = ifelse(is.null(egger_res), NA, egger_res$intercept),
    se_egger = ifelse(is.null(egger_res), NA, egger_res$intercept_se),
    pval_egger = ifelse(is.null(egger_res), NA, egger_res$intercept_pval)
  )

  # Pleiotropy
  pleio_res <- tryCatch(mr_pleiotropy_test(harmonized_dat), error = function(e) NULL)
  if (!is.null(pleio_res) && nrow(pleio_res) > 0) {
    pleio_res$trait <- lipid_name
    pleio_res$method <- "MR Egger"
    pleio_results_list[[lipid_name]] <- pleio_res
  }

  # Heterogeneity
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
          presso_global_pval = presso_res$`MR-PRESSO results`$`Global Test`$`P-value`,
          distortion_pval = ifelse(!is.null(presso_res$`Distortion Test`$`P-value Distortion`),
                                   presso_res$`Distortion Test`$`P-value Distortion`, NA),
          outlier_count = length(presso_res$`MR-PRESSO results`$`Outliers Indices`)
        )
      }
    }, error = function(e) {
      cat("MR-PRESSO failed for:", lipid_name, ":", conditionMessage(e), "\n")
    })
  }

  cat("Completed:", lipid_name, "\n")
}

# Combine all results
combined_results <- bind_rows(mr_results_list)
combined_egger <- bind_rows(egger_intercepts_list)
combined_pleio <- bind_rows(pleio_results_list)
combined_hetero <- bind_rows(hetero_results_list)
combined_presso <- bind_rows(presso_results_list)

# Merge all components
final_combined <- combined_results
if (nrow(combined_egger) > 0) {
  final_combined <- left_join(final_combined, combined_egger, by = "trait")
}
if (nrow(combined_pleio) > 0) {
  pleio_df <- combined_pleio %>%
    select(trait, method, egger_intercept, se, pval) %>%
    rename(pleio_egger_intercept = egger_intercept, pleio_se = se, pleio_pval = pval)
  final_combined <- left_join(final_combined, pleio_df, by = c("trait", "method"))
}
if (nrow(combined_hetero) > 0) {
  final_combined <- left_join(final_combined, combined_hetero %>% select(trait, method, Q, Q_df, Q_pval), by = c("trait", "method"))
}
if (nrow(combined_presso) > 0) {
  final_combined <- left_join(final_combined, combined_presso, by = "trait")
}

# Save results
write.csv(final_combined, "final_combined_reverse_mr_results_CAD_as_exposure_with_presso.csv", row.names = FALSE)

cat("Reverse MR with MR-PRESSO complete!\n")
