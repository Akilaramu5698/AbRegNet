library(sva)

expr_all <- read.csv("all_GEO_merged_logCPM.csv", row.names = 1, check.names = FALSE)
expr_all <- as.matrix(expr_all); mode(expr_all) <- "numeric"

# --- Global ComBat run: dataset-wise batch, 90 samples across 7 GEO batches ---
meta <- read.csv("metadata_global_dataset_wise.csv", check.names = FALSE)
expr <- expr_all[, meta$Run]
stopifnot(all(colnames(expr) == meta$Run))

batch <- factor(meta$`GEO ID`)
mod <- model.matrix(~ Antibiotic, data = meta)   # preserves biological signal

global_combat <- ComBat(dat = expr, batch = batch, mod = mod, par.prior = TRUE)
write.csv(global_combat, "global_ComBat_logCPM.csv")

# --- Tigecycline: single batch, no correction possible/needed ---
tige_meta <- read.csv("metadata_tigecycline_only.csv", check.names = FALSE)
tige_expr <- expr_all[, tige_meta$Run]
write.csv(tige_expr, "tigecycline_logCPM_no_ComBat.csv")

# --- Downstream: slice the globally-corrected matrix per antibiotic for each GRN ---
for (drug in c("Amikacin","Ciprofloxacin","Meropenem","Colistin/Polymyxin")) {
  samples <- meta$Run[meta$Antibiotic %in% c(drug, "Control")]
  grn_expr <- global_combat[, samples]
  write.csv(grn_expr, paste0(gsub("/","-",drug), "_for_GRN.csv"))
}


# --- Merge global ComBat matrix + tigecycline (uncorrected) into one 96-sample file ---

global_combat <- read.csv("global_ComBat_logCPM.csv", row.names = 1, check.names = FALSE)
tige_expr     <- read.csv("tigecycline_logCPM_no_ComBat.csv", row.names = 1, check.names = FALSE)

# both share the same gene rows (same row.names / gene IDs), so cbind on columns
stopifnot(all(rownames(global_combat) == rownames(tige_expr)))

merged_96 <- cbind(global_combat, tige_expr)

cat("Merged matrix:", nrow(merged_96), "genes x", ncol(merged_96), "samples\n")

write.csv(merged_96, "all96_merged_logCPM_ComBat.csv")
