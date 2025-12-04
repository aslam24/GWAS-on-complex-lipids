# Load packages
library(data.table)
library(coloc)

# Step 1: Load T2D dataset
T2D <- fread("T2D_Xue_et_al_2018.txt")
T2D <- T2D[!duplicated(SNP)]
T2D <- T2D[complete.cases(T2D[, .(b, se, frq_A1, N)])]  # Remove any NA rows

# Step 2: Load lipid datasets
lipid_files <- list.files("lipid_files_ZPR1", full.names = TRUE, pattern = "\\.txt$")
lipid_list <- lapply(lipid_files, fread)

# Step 3: Prepare to store coloc results
coloc_results <- vector("list", length(lipid_list))

# Step 4: Loop through each lipid trait
for (i in seq_along(lipid_list)) {
  lipid_data <- lipid_list[[i]]
  lipid_data <- lipid_data[!duplicated(SNP)]
  lipid_data <- lipid_data[complete.cases(lipid_data[, .(SNP, ES, SE, MAF, N)])]
  
  # Merge datasets on SNP
  merged <- merge(T2D, lipid_data, by = "SNP")
  
  # Filter missing values
  merged <- merged[complete.cases(merged[, .(b, se, ES, SE, MAF, frq_A1)])]
  
  # Skip if too few SNPs
  if (nrow(merged) < 100) {
    message(sprintf("Skipping trait %d due to low SNP overlap", i))
    next
  }
  
  # Define datasets
  t2d_dataset <- list(
    beta = merged$b,
    varbeta = merged$se^2,
    snp = merged$SNP,
    MAF = merged$frq_A1,
    N = as.integer(median(merged$N.x)),  # robust estimate
    type = "cc",
    s = 0.1  # case proportion (adjust if known)
  )
  
  lipid_dataset <- list(
    beta = merged$ES,
    varbeta = merged$SE^2,
    snp = merged$SNP,
    MAF = merged$MAF,
    N = merged$N.y[1],
    type = "quant"
  )
  
  # Run coloc
  coloc_results[[i]] <- coloc.abf(dataset1 = t2d_dataset, dataset2 = lipid_dataset)
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

# Step 6: Top results by posterior probability of colocalization (PP.H4)
top_hits <- summary_df[order(-PP.H4)]
head(top_hits, 10)
