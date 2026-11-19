setwd("F:/current-work/2nd.obj/normalization")

library(edgeR)

# Find all GEO count CSV files
files <- list.files(
  pattern = "_counts\\.csv$",
  full.names = TRUE
)

for (file in files) {
  
  cat("\nProcessing:", basename(file), "\n")
  
  # Read CSV
  merged <- read.csv(
    file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  # Check NA
  cat("NA values:", sum(is.na(merged)), "\n")
  

  
  # Convert to numeric
  merged <- as.data.frame(merged)
  merged[] <- lapply(merged, as.numeric)
  
  # Create DGEList
  dge <- DGEList(counts = merged)
  
  # Upper-quartile normalization
  dge <- calcNormFactors(
    dge,
    method = "upperquartile"
  )
  
  # Log2 CPM
  logCPM <- cpm(
    dge,
    log = TRUE,
    prior.count = 1
  )
  
  # Create output filename
  output_file <- sub(
    "_counts\\.csv$",
    "_upperquartile_logCPM.csv",
    file
  )
  
  # Save
  write.csv(
    logCPM,
    output_file
  )
  
  cat("Saved:", basename(output_file), "\n")
}



setwd("F:/current-work/2nd.obj/normalization")

library(dplyr)

files <- list.files(
  pattern = "_upperquartile_logCPM\\.csv$",
  full.names = TRUE
)

# Read first file
merged_logCPM <- read.csv(
  files[1],
  row.names = 1,
  check.names = FALSE
)

# Add remaining files column-wise
for (file in files[-1]) {
  
  temp <- read.csv(
    file,
    row.names = 1,
    check.names = FALSE
  )
  
  # Check that genes are the same
  if (!identical(rownames(merged_logCPM), rownames(temp))) {
    stop(
      paste(
        "Gene mismatch in:",
        basename(file)
      )
    )
  }
  
  merged_logCPM <- cbind(
    merged_logCPM,
    temp
  )
}

# Save merged file
write.csv(
  merged_logCPM,
  "all_GEO_merged_logCPM.csv"
)

cat("Merged", length(files), "GEO files\n")
cat("Genes:", nrow(merged_logCPM), "\n")
cat("Samples:", ncol(merged_logCPM), "\n")


metadata <- read.table(
  "meta.data.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)


