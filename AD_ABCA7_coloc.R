library(data.table)
library(coloc)

# Step 1: Load AD dataset
AD <- fread("AD_2022_Grch37.txt", header = TRUE)

# Remove duplicates and incomplete rows
AD <- AD[!duplicated(variant_id)]
AD <- AD[complete.cases(AD[, .(beta, standard_error, effect_allele_frequency, n_cases, n_controls)])]

# Rename columns to match coloc expected input
setnames(AD, c("variant_id", "beta", "standard_error", "effect_allele_frequency", "n_cases", "n_controls"),
            c("SNP", "b", "se", "frq_A1", "cases", "controls"))

# Calculate total sample size and case proportion
AD[, N := cases + controls]
AD[, s := cases / N]

# Step 2: Load lipid datasets
lipid_files <- list.files("lipid_list_ABCA7", full.names = TRUE, pattern = "\\.txt$")
lipid_list <- lapply(lipid_files, fread)

# Step 3: Prepare to store coloc results
coloc_results <- vector("list", length(lipid_list))

# Step 4: Loop through each lipid trait
for (i in seq_along(lipid_list)) {
  lipid_data <- lipid_list[[i]]
  lipid_data <- lipid_data[!duplicated(SNP)]
  lipid_data <- lipid_data[complete.cases(lipid_data[, .(SNP, ES, SE, MAF, N)])]
  
  # Merge datasets on SNP
  merged <- merge(AD, lipid_data, by = "SNP")
  
  # Filter missing values
  merged <- merged[complete.cases(merged[, .(b, se, ES, SE, MAF, frq_A1)])]
  
  # Skip if too few SNPs
  if (nrow(merged) < 100) {
    message(sprintf("Skipping trait %d due to low SNP overlap", i))
    next
  }
  
  # Define AD dataset for coloc (case-control)
  AD_dataset <- list(
    beta = merged$b,
    varbeta = merged$se^2,
    snp = merged$SNP,
    MAF = merged$frq_A1,
    N = as.integer(median(merged$N)),   # median total sample size
    type = "cc",
    s = median(merged$s)                # case fraction
  )
  
  # Define lipid dataset (quantitative)
  lipid_dataset <- list(
    beta = merged$ES,
    varbeta = merged$SE^2,
    snp = merged$SNP,
    MAF = merged$MAF,
    N = merged$N.y[1],
    type = "quant"
  )
  
  # Run coloc
  coloc_results[[i]] <- coloc.abf(dataset1 = AD_dataset, dataset2 = lipid_dataset)
}

# Step 5: Summarize results
summary_df <- rbindlist(
  lapply(seq_along(coloc_results), function(i) {
    res <- coloc_results[[i]]
    trait_name <- if (!is.null(lipid_list[[i]])) unique(lipid_list[[i]]$Lipid) else NA
    if (!is.null(res) && !is.null(res$summary)) {
      data.table(
        index = i,
        trait = trait_name,
        nsnps = res$summary["nsnps"],
        PP.H0 = res$summary["PP.H0.abf"],
        PP.H1 = res$summary["PP.H1.abf"],
        PP.H2 = res$summary["PP.H2.abf"],
        PP.H3 = res$summary["PP.H3.abf"],
        PP.H4 = res$summary["PP.H4.abf"]
      )
    } else {
      data.table(
        index = i,
        trait = trait_name,
        nsnps = NA, PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA
      )
    }
  }),
  fill = TRUE
)

# Step 6: Top results by PP.H4
top_hits <- summary_df[order(-PP.H4)]
head(top_hits, 10)

# Step 7: Save results
fwrite(summary_df, file = "coloc_AD_ABCA7_summary.tsv", sep = "\t")
fwrite(head(top_hits, 10), file = "coloc_AD_ABCA7_top_hits.tsv", sep = "\t")
