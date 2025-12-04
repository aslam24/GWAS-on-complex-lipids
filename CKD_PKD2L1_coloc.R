library(data.table)
library(coloc)

###------------------------------------------------------------
### Step 1: Load CKD dataset
###------------------------------------------------------------

CKD <- fread("CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt")

# Remove duplicates and incomplete rows
CKD <- CKD[!duplicated(RSID)]
CKD <- CKD[complete.cases(CKD[, .(Effect, StdErr, Freq1, n_total_sum)])]

# Rename columns
setnames(
  CKD,
  c("RSID", "Effect", "StdErr", "Freq1", "n_total_sum"),
  c("SNP", "b", "se", "frq_A1", "N")
)

###------------------------------------------------------------
### Step 2: Load all lipid datasets
###------------------------------------------------------------

lipid_files <- list.files("lipid_list_PKD2L1", full.names = TRUE, pattern = "\\.txt$")
lipid_list  <- lapply(lipid_files, fread)
names(lipid_list) <- gsub(".txt$", "", basename(lipid_files))

###------------------------------------------------------------
### Step 3: Store coloc results
###------------------------------------------------------------

coloc_results <- vector("list", length(lipid_list))

###------------------------------------------------------------
### Step 4: Loop through lipid traits and run coloc
###------------------------------------------------------------

for (i in seq_along(lipid_list)) {

  lipid_data <- lipid_list[[i]]

  lipid_data <- lipid_data[!duplicated(SNP)]
  lipid_data <- lipid_data[complete.cases(lipid_data[, .(SNP, ES, SE, MAF, N)])]

  # Merge datasets
  merged <- merge(CKD, lipid_data, by = "SNP")
  merged <- merged[complete.cases(merged[, .(b, se, ES, SE, MAF, frq_A1)])]

  # Skip insufficient SNP overlap
  if (nrow(merged) < 100) {
    message(sprintf("Skipping %s (index %d): too few SNPs (%d)", names(lipid_list)[i], i, nrow(merged)))
    next
  }

  # Define dataset 1 (CKD)
  CKD_dataset <- list(
    beta = merged$b,
    varbeta = merged$se^2,
    snp = merged$SNP,
    MAF = merged$frq_A1,
    N = as.integer(median(merged$N)),
    type = "quant"
  )

  # Dataset 2 (lipid)
  lipid_dataset <- list(
    beta = merged$ES,
    varbeta = merged$SE^2,
    snp = merged$SNP,
    MAF = merged$MAF,
    N = merged$N.y[1],
    type = "quant"
  )

  # Run coloc
  coloc_results[[i]] <- coloc.abf(dataset1 = CKD_dataset, dataset2 = lipid_dataset)
}

###------------------------------------------------------------
### Step 5: Summarize coloc results
###------------------------------------------------------------

summary_df <- rbindlist(
  lapply(seq_along(coloc_results), function(i) {
    res <- coloc_results[[i]]

    trait_name <- names(lipid_list)[i]

    if (!is.null(res) && !is.null(res$summary)) {
      return(data.table(
        trait = trait_name,
        nsnps = as.numeric(res$summary["nsnps"]),
        PP.H0 = as.numeric(res$summary["PP.H0.abf"]),
        PP.H1 = as.numeric(res$summary["PP.H1.abf"]),
        PP.H2 = as.numeric(res$summary["PP.H2.abf"]),
        PP.H3 = as.numeric(res$summary["PP.H3.abf"]),
        PP.H4 = as.numeric(res$summary["PP.H4.abf"])
      ))
    } else {
      return(data.table(
        trait = trait_name,
        nsnps = NA, PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA
      ))
    }
  }),
  fill = TRUE
)

###------------------------------------------------------------
### Step 6: Sort by strongest colocalization (PP.H4)
###------------------------------------------------------------

top_hits <- summary_df[order(-PP.H4)][1:10]

###------------------------------------------------------------
### Step 7: Save results
###------------------------------------------------------------

fwrite(summary_df, "coloc_CKD_PKD2L1_summary.tsv", sep = "\t")
fwrite(top_hits, "coloc_CKD_PKD2L1_top_hits.tsv", sep = "\t")
