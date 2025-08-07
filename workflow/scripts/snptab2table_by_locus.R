# Load arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
files <- args[1:(length(args) - 1)]
locus_tags <- unlist(strsplit(args[length(args)], " "))

# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(forcats)
library(openxlsx)

# Read SNP tables
snp_list <- lapply(files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  df$Sample <- gsub(".snps.tab", "", basename(f))
  df
})
all_snps <- bind_rows(snp_list)

# Filter by Locus Tag
filtered <- all_snps %>% filter(LOCUS_TAG %in% locus_tags)
cat("Number of SNPs in selected loci:", nrow(filtered), "\n")

# Generate LOCUS_POS column
filtered <- filtered %>% mutate(LOCUS_POS = paste0(LOCUS_TAG, "_", POS))

# Simplify effect categories
filtered$EFFECT_SIMPLE <- case_when(
  grepl("missense_variant", filtered$EFFECT) ~ "missense_variant",
  grepl("synonymous_variant", filtered$EFFECT) ~ "synonymous_variant",
  grepl("frameshift", filtered$EFFECT) ~ "frameshift_variant",
  grepl("stop", filtered$EFFECT) ~ "stop_variant",
  grepl("splice", filtered$EFFECT) ~ "splice_variant",
  TRUE ~ "other"
)

# Create presence matrix
mat <- filtered %>%
  mutate(value = 1) %>%
  select(Sample, LOCUS_POS, value) %>%
  pivot_wider(names_from = LOCUS_POS, values_from = value, values_fill = 0)

mat_long <- mat %>%
  pivot_longer(cols = -Sample, names_to = "LOCUS_POS", values_to = "Present")

# Merge annotations
annotations <- filtered %>%
  select(LOCUS_POS, LOCUS_TAG, GENE, TYPE, EFFECT_SIMPLE) %>%
  distinct()

mat_long <- mat_long %>%
  left_join(annotations, by = "LOCUS_POS") %>%
  mutate(
    TYPE = ifelse(is.na(TYPE), "none", TYPE),
    EFFECT_SIMPLE = ifelse(is.na(EFFECT_SIMPLE), "none", EFFECT_SIMPLE)
  )

# Load metadata
strain_info <- read.xlsx("output/variant_calling/snptab_viewer/StrainList.xlsx")
name_map <- read.xlsx("output/variant_calling/snptab_viewer/name_map.xlsx") %>%
  rename(SRR_ID = old_name, patient = new_name)

# Map SRR IDs to patient names, keep EV_* samples as is
mat_long <- mat_long %>%
  left_join(name_map, by = c("Sample" = "SRR_ID")) %>%
  mutate(
    patient_final = ifelse(!is.na(patient), patient, Sample)  # SRR mapped to patient, EV stays
  )

# Join full metadata
mat_long <- left_join(mat_long, strain_info, by = c("patient_final" = "patient"))

# Final table
output_table <- mat_long %>%
  select(Sample, LOCUS_POS, LOCUS_TAG, GENE, Present, TYPE, EFFECT_SIMPLE, cohort, bacteremiaduration) %>%
  arrange(Sample, LOCUS_POS)

# Save
output_dir <- "output/variant_calling/snptab_viewer"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv(output_table, file.path(output_dir, "snp_by_locus_data.tsv"))
write.xlsx(output_table, file.path(output_dir, "snp_by_locus_data.xlsx"))