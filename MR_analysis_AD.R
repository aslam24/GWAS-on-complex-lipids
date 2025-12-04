library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)
library(RadialMR)
library(ieugwasr)

plink_bin <- genetics.binaRies::get_plink_binary()
bfile_ref <- "EUR"  # Update as needed

exposure_files <- list.files("AD_MR", pattern = "\\.txt$", full.names = TRUE)
AD <- fread("AD_2022_Grch37.txt")

mr_results_list <- list()
mr_results_excl_list <- list()  # for radial MR filtered results
pleio_results_list <- list()
hetero_results_list <- list()
egger_intercepts_list <- list()

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
  
  AD_sub <- AD[AD$variant_id %in% exposure_dat$SNP, ]
  if (nrow(AD_sub) == 0) {
    cat("No overlapping SNPs in AD outcome. Skipping.\n")
    next
  }
  
  outcome_dat <- format_data(
    as.data.frame(AD_sub),
    type = "outcome",
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
  
  harmonized_dat$outcome <- "AD"
  
  # Run MR on full harmonized data
  mr_res <- tryCatch({
    mr(harmonized_dat)
  }, error = function(e) {
    cat("MR failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(mr_res)) next
  
  # Add trait and nsnp info
  mr_res$trait <- lipid_name
  mr_res$nsnp <- nrow(harmonized_dat)
  
  # Collect pleio and heterogeneity test results for original MR only
  pleio_res <- tryCatch(mr_pleiotropy_test(harmonized_dat), error = function(e) NULL)
  hetero_res <- tryCatch(mr_heterogeneity(harmonized_dat), error = function(e) NULL)
  egger_res <- tryCatch(mr_egger_regression(harmonized_dat), error = function(e) NULL)
  
  mr_results_list[[lipid_name]] <- mr_res
  
  egger_intercepts_list[[lipid_name]] <- data.frame(
    trait = lipid_name,
    intercept_egger = ifelse(is.null(egger_res), NA, egger_res$intercept),
    se_egger = ifelse(is.null(egger_res), NA, egger_res$intercept_se),
    pval_egger = ifelse(is.null(egger_res), NA, egger_res$intercept_pval)
  )
  
  if (!is.null(pleio_res)) {
    pleio_res$trait <- lipid_name
    pleio_res$method <- "MR Egger"
    pleio_results_list[[lipid_name]] <- pleio_res
  }
  
  if (!is.null(hetero_res)) {
    hetero_res$trait <- lipid_name
    hetero_results_list[[lipid_name]] <- hetero_res
  }
  
  # --- Radial MR outlier detection and filtering if heterogeneity significant ---
  if (!is.null(hetero_res) && any(hetero_res$Q_pval < 0.05, na.rm = TRUE)) {
    cat("Significant heterogeneity detected for", lipid_name, "- running Radial MR outlier detection\n")
    
    if (nrow(harmonized_dat) >= 3) {
      for_rad <- format_radial(
        BXG = harmonized_dat$beta.exposure,
        BYG = harmonized_dat$beta.outcome,
        seBXG = harmonized_dat$se.exposure,
        seBYG = harmonized_dat$se.outcome,
        RSID = harmonized_dat$SNP
      )
      
      ivw_radial_res <- tryCatch(ivw_radial(for_rad), error = function(e) {
        cat("ivw_radial error:", conditionMessage(e), "\n")
        return(NULL)
      })
      
      egger_radial_res <- tryCatch(egger_radial(for_rad), error = function(e) {
        cat("egger_radial error:", conditionMessage(e), "\n")
        return(NULL)
      })
      
      # Safe check for outliers extraction
      if (!is.null(ivw_radial_res)) {
        if (is.list(ivw_radial_res) && "outliers" %in% names(ivw_radial_res)) {
          if (!is.null(ivw_radial_res$outliers) && nrow(ivw_radial_res$outliers) > 0) {
            outliers <- ivw_radial_res$outliers$SNP
            cat("Removing", length(outliers), "outliers from harmonized data\n")
            harmonized_no_outliers <- harmonized_dat[!harmonized_dat$SNP %in% outliers, ]
            
            if (nrow(harmonized_no_outliers) >= 3) {
              mr_res_excl <- tryCatch({
                mr(harmonized_no_outliers)
              }, error = function(e) {
                cat("MR on filtered data failed:", conditionMessage(e), "\n")
                return(NULL)
              })
              
              # Collect pleiotropy and heterogeneity tests on filtered data
              pleio_res_excl <- tryCatch(mr_pleiotropy_test(harmonized_no_outliers), error = function(e) NULL)
              hetero_res_excl <- tryCatch(mr_heterogeneity(harmonized_no_outliers), error = function(e) NULL)
              
              if (!is.null(mr_res_excl)) {
                if (!is.null(pleio_res_excl)) {
                  pleio_res_excl_renamed <- pleio_res_excl %>%
                    rename(
                      pleio_egger_intercept_filtered = egger_intercept,
                      pleio_se_filtered = se,
                      pleio_pval_filtered = pval
                    ) %>%
                    select(method, pleio_egger_intercept_filtered, pleio_se_filtered, pleio_pval_filtered)
                } else {
                  pleio_res_excl_renamed <- NULL
                }
                
                if (!is.null(hetero_res_excl)) {
                  hetero_res_excl_renamed <- hetero_res_excl %>%
                    rename(
                      Q_filtered = Q,
                      Q_df_filtered = Q_df,
                      Q_pval_filtered = Q_pval
                    ) %>%
                    select(method, Q_filtered, Q_df_filtered, Q_pval_filtered)
                } else {
                  hetero_res_excl_renamed <- NULL
                }
                
                mr_res_excl <- mr_res_excl %>%
                  { if (!is.null(pleio_res_excl_renamed)) left_join(., pleio_res_excl_renamed, by = "method") else . } %>%
                  { if (!is.null(hetero_res_excl_renamed)) left_join(., hetero_res_excl_renamed, by = "method") else . }
                
                if (!"trait" %in% colnames(mr_res_excl)) mr_res_excl$trait <- paste0(lipid_name, "_excl_outliers")
                if (!"nsnp" %in% colnames(mr_res_excl)) mr_res_excl$nsnp <- nrow(harmonized_no_outliers)
                
                mr_results_excl_list[[paste0(lipid_name, "_excl_outliers")]] <- mr_res_excl
              }
            } else {
              cat("Not enough SNPs left after outlier removal for", lipid_name, "\n")
            }
          } else {
            cat("No outliers detected by Radial MR for", lipid_name, "\n")
          }
        } else {
          cat("No outliers detected by Radial MR for", lipid_name, "\n")
        }
      } else {
        cat("ivw_radial returned NULL for", lipid_name, "\n")
      }
    } else {
      cat("Not enough SNPs (less than 3) to run Radial MR for", lipid_name, "\n")
    }
  }
  
  cat("Completed:", lipid_name, "\n")
}

# Combine all MR results
combined_results <- do.call(rbind, mr_results_list)
combined_results_excl <- if(length(mr_results_excl_list) > 0) do.call(rbind, mr_results_excl_list) else NULL
combined_pleio <- do.call(rbind, pleio_results_list)
combined_hetero <- do.call(rbind, hetero_results_list)
combined_egger <- do.call(rbind, egger_intercepts_list)

library(dplyr)

# Rename pleio columns for merging
pleio_df <- combined_pleio %>%
  select(trait, method, egger_intercept, se, pval) %>%
  rename(
    pleio_egger_intercept = egger_intercept,
    pleio_se = se,
    pleio_pval = pval
  )

# Merge combined_results with pleio results
merged1 <- combined_results %>%
  left_join(pleio_df, by = c("trait", "method"))

# Merge with heterogeneity results
final_combined <- merged1 %>%
  left_join(
    combined_hetero %>%
      select(trait, method, Q, Q_df, Q_pval),
    by = c("trait", "method")
  )

# Append filtered MR results with additional hetero and pleio cols (radial MR outlier exclusion)
if (!is.null(combined_results_excl) && nrow(combined_results_excl) > 0) {
  final_combined <- bind_rows(final_combined, combined_results_excl)
}

# View combined data
head(final_combined)

# Export combined table to file
fwrite(final_combined, "MR_results_full_combined.txt", sep = "\t")
