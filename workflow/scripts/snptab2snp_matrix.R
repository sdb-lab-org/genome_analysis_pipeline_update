#!/usr/bin/env Rscript

# ---- Load arguments ----
args <- commandArgs(trailingOnly = TRUE)
files <- args[1:(length(args)-2)]
genes_arg <- args[length(args)-1]
output_file <- args[length(args)]

genes_of_interest <- unlist(strsplit(genes_arg, " "))

cat("Genes of interest:", genes_of_interest, "\n")
cat("Output file:", output_file, "\n")

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ---- Read all SNP tables ----
all_snps <- lapply(files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  df$Sample <- gsub(".snps.tab$", "", basename(f))
  df
})

# Standardize POS column to numeric
all_snps <- lapply(all_snps, function(df) {
  df$POS <- as.numeric(df$POS)
  df
})

# Bind all together
all_snps_df <- bind_rows(all_snps)

# ---- Filter by genes of interest ----
filtered <- all_snps_df %>% filter(GENE %in% genes_of_interest)

if (nrow(filtered) == 0) {
  cat("No SNPs found for the selected genes.\n")
  # Write an empty table with headers
  write_tsv(data.frame(GENE_POS=character(), GENE=character(), EFFECT=character()), output_file)
  quit(save="no", status=0)
}

# Create a unique key for each SNP
filtered <- filtered %>% mutate(GENE_POS = paste0(GENE, "_", POS))

# Pivot to presence/absence
presence_df <- filtered %>%
  mutate(Present = 1) %>%
  select(Sample, GENE_POS, GENE, EFFECT, Present) %>%
  pivot_wider(
    id_cols = c(GENE_POS, GENE, EFFECT),
    names_from = Sample,
    values_from = Present,
    values_fill = 0
  )

# Ensure every sample appears as a column
all_samples <- unique(basename(gsub(".snps.tab$", "", files)))
for (s in all_samples) {
  if (!(s %in% colnames(presence_df))) {
    presence_df[[s]] <- 0
  }
}

# Reorder columns: GENE_POS, GENE, EFFECT, then samples in original order
presence_df <- presence_df %>%
  select(GENE_POS, GENE, EFFECT, all_of(all_samples))

# ---- Write output ----
write_tsv(presence_df, output_file)

cat("Matrix written to", output_file, "\n")