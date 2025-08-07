#############################################################################################
##                        Snakemake pipeline for Genome analysis                           ##
#############################################################################################

########################################
## To run the pipeline use this command#
########################################

'''
snakemake -s workflow/Snakefile.py --workflow-profile workflow/profiles/ga_pipeline/ -n
snakemake -s workflow/Snakefile.py --profile workflow/profiles/default -n
'''

########################################
## Python packages                    ##
########################################

import pandas as pd
import glob
import os
import yaml

########################################
## Configuration                      ##
########################################

configfile: "config/config.yaml" 

WORKING_DIR         = config["working_dir"]
OUTPUT_DIR          = config["output_dir"]
DATABASE_DIR        = config["database_dir"]
CGMLST_SCHEMA_DIR   = config["cgmlst_schema"] 
PROJECT_DIR         = config["project_dir"]
PROJECT_NAME        = config["project_name"]
ASSEMBLER           = config["assembler"]
ANNOTATOR           = config["annotator"]

########################################
## Samples                            ##
########################################

## Raw read data
SAMPLES         = config["samples"]

## External genomes 
GENOMES         = config["genomes"]

ALL_GENOMES     = sorted(set(SAMPLES + GENOMES))

# Reference
REF_SAMPLE      = config["refg"][0]

GWAS_GENOMES = [g for g in ALL_GENOMES if g != REF_SAMPLE]

print("SAMPLES:", SAMPLES) 
print("GENOMES:", GENOMES)
print("ALL_GENOMES:", ALL_GENOMES)
print("REF_SAMPLE:", REF_SAMPLE)
print("GWAS_GENOMES:", GWAS_GENOMES)

########################################
## Rules                              ##
########################################

# QC
include: "rules/fastp.smk"
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"

## Assembly
if ASSEMBLER == "spades":
    include: "rules/spades.smk"
elif ASSEMBLER == "unicycler":
    include: "rules/unicycler.smk"
else:
    raise ValueError(f"Unknown assembler specified: {ASSEMBLER}")
include: "rules/quast.smk"

# Copy external genomes
include: "rules/copy_genomes.smk"

## Annotation
if ANNOTATOR == "prokka":
    include: "rules/prokka.smk"
elif ANNOTATOR == "bakta":
    include: "rules/bakta.smk"
else:
    raise ValueError(f"Unknown aannotator specified: {ANNOTATOR}")
include: "rules/emapper_kegganog.smk"


# QC summary
include: "rules/summarize_genomes.smk" 

# Pangenome analysis
include: "rules/pirate.smk"
include: "rules/anvio.smk"

# Typing
include: "rules/mlst_spatyper.smk"
include: "rules/chewbacca.smk" 

# AMR/virulence genes screening
include: "rules/abricate.smk"

# SNP analysis 
include: "rules/refg_gbk.smk"
include: "rules/snippy.smk"

# GWAS analysis 
include: "rules/gwas.smk"

# Mapping
include: "rules/align_reads_to_reference.smk"

# Filter genes of interest
include: "rules/genes_of_interest.smk"

########################################
## Desired outputs                    ##
########################################

# Assembly and QC 
QC_MULTIQC                      = [OUTPUT_DIR + "01_qc/multiqc/multiqc_report.html"]
QC_FASTP                        = expand(OUTPUT_DIR + "01_qc/fastp/{sample}_trimmed_R1.fastq.gz", sample=SAMPLES)
ASSEMBLY_CONTIGS                = expand(OUTPUT_DIR + "08_temp/temp_fasta/{sample}.fasta", sample=SAMPLES)
QC_QUAST                        = expand(OUTPUT_DIR + "01_qc/quast/{sample}/report.txt", sample=SAMPLES)

# Copy external genomes
COPY_EXTERNAL_GENOMES           = expand(OUTPUT_DIR + "08_temp/temp_fasta/{genome}.fasta", genome=GENOMES)

# Annotation 
ANNOTATION_COPY_GFF      = expand(OUTPUT_DIR + "08_temp/temp_gff/{all_genomes}.gff", all_genomes=ALL_GENOMES)

# QC summary 
QC_SUMMARY                      = [OUTPUT_DIR + "01_qc/summary_qc.tsv"] 

# Annotation KEGGANOG/EMAPPER
ANNOTATION_EMAPPER              = expand(OUTPUT_DIR + "{all_genomes}/annotation/emapper/", all_genomes=ALL_GENOMES)
KEGGANOG                        = [OUTPUT_DIR + "02_kegganog/heatmap_figure.png"]
ANNOTATION_COPY_DEC_GFF         = expand(OUTPUT_DIR + "08_temp/temp_decorated_gff/{all_genomes}.emapper.decorated.gff", all_genomes=ALL_GENOMES)


# Pangenome analysis PIRATE 
PANGENOME_PIRATE                = [OUTPUT_DIR + "03_pangenome/pirate/core_alignment.fasta"]
PANGENOME_TREE                  = [OUTPUT_DIR + "03_pangenome/pirate/core_alignment.fasta.treefile"]
PANGENOME_CORE_GENES_DISTANCE   = [OUTPUT_DIR + "03_pangenome/pirate/core_gene_snp_distances.tsv"]

# Typing
TYPING_MLST                     = [OUTPUT_DIR + "06_typing/mlst/mlst.tsv"]
TYPING_SPATYPER                 = expand(OUTPUT_DIR + "06_typing/spaTyper/{all_genomes}_spatype.txt", all_genomes=ALL_GENOMES)
CGMLST_PREP_EXTERNAL_DB         = [CGMLST_SCHEMA_DIR + "PreparedSchema/schema.json"] 
CGMLST_ALLELE_CALL              = [OUTPUT_DIR + "06_typing/chewbbaca/results_allele_call.json"]

# AMR/virulence genes screening 
AMR_NCBI                        = expand(OUTPUT_DIR + "07_amr/{all_genomes}.ncbi.tsv", all_genomes=ALL_GENOMES) 
AMR_NCBI_SUMMARY                = [OUTPUT_DIR + "07_amr/amr_ncbi_summary.tab"]


# Pangenome analysis ANVI'O 
ANVIO_FASTA_REFORMATTED         = expand(OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.simplified.fasta", all_genomes=ALL_GENOMES) 
ANVIO_GENE_CALLS                = expand(OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}_gene_calls.txt", all_genomes=ALL_GENOMES) 
ANVIO_ANNOTATION                = expand(OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}_gene_annot.txt", all_genomes=ALL_GENOMES)
ANVIO_DB                        = expand(OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}.db", all_genomes=ALL_GENOMES) 
ANVIO_DONE_DB                   = expand(OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}.annotation.done", all_genomes=ALL_GENOMES) 
ANVIO_EXTERNAL_GENOMES_FILE     = [OUTPUT_DIR + "03_pangenome/anvio/external-genomes.txt"]
ANVIO_STORAGE_DB                = [OUTPUT_DIR + "03_pangenome/anvio/anvio_storage-GENOMES.db"]
ANVIO_PANGENOME_DB              = [OUTPUT_DIR + "03_pangenome/anvio/anvio_Pangenome-PAN.txt"]
ANVIO_ANI                       = [OUTPUT_DIR + "03_pangenome/anvio/ANI/genome_similarity.txt"]


# Variant calling SNIPPY 
VC_REFG_GBK                     = [WORKING_DIR + "reference/refg.gbk"]
VC_SNIPPY                       = expand(OUTPUT_DIR + "04_variant_calling/snippy/{all_genomes}/snps.vcf", all_genomes=ALL_GENOMES)
VC_COPY_VCF                     = expand(OUTPUT_DIR + "08_temp/temp_vcf/{all_genomes}.snps.vcf", all_genomes=ALL_GENOMES)
VC_COPY_VCF_GZ                  = expand(OUTPUT_DIR + "08_temp/temp_vcf_gz/{all_genomes}.snps.vcf.gz", all_genomes=ALL_GENOMES) 
VC_COPY_TAB                     = expand(OUTPUT_DIR + "08_temp/temp_tab/{all_genomes}.snps.tab", all_genomes=ALL_GENOMES)
VC_SNIPPY_CORE                  = [OUTPUT_DIR + "04_variant_calling/snippy-core/core.full.aln"]
VC_CORE_TREE                    = [OUTPUT_DIR + "04_variant_calling/snippy-core/iqtree.log"]
VC_VCF_HEATMAP                  = [OUTPUT_DIR + "04_variant_calling/vcf_viewer/heatmap_output.html"]
VC_SNP_DISTANCE                 = [OUTPUT_DIR + "04_variant_calling/snippy-core/snp_distance_matrix.tsv"]

# GWAS 
CSI                             = expand(OUTPUT_DIR + "08_temp/temp_vcf_gz/{gwas_genomes}.snps.vcf.gz.csi", gwas_genomes=GWAS_GENOMES) 
VCF_MERGED                      = [OUTPUT_DIR + "05_gwas/merged.vcf.gz"]
PHYLO_DIST                      = [OUTPUT_DIR + "05_gwas/phylogeny_dists.tsv"]


# Mapping
BAM                             = expand(OUTPUT_DIR + "{sample}/mapping/{sample}.sorted.bam", sample=SAMPLES)

# Filter genes of interests
PIRATE_GENE_TABLE               = [OUTPUT_DIR + "03_pangenome/genes_of_interest_presence_absence.tsv"]
SNP_GENES_OF_INTEREST           = [OUTPUT_DIR + "04_variant_calling/genes_of_interest_snp_matrix.tsv"]

########################################
## Pipeline                           ##
########################################

rule all:
    input:
        QC_MULTIQC,
        ASSEMBLY_CONTIGS,
        QC_QUAST,
        COPY_EXTERNAL_GENOMES,
        ANNOTATION_COPY_GFF, 
        ##### Optional steps depending on config emapper_kegganog = true
        *(ANNOTATION_EMAPPER if config["run"].get("emapper_kegganog") else []),
        *(ANNOTATION_COPY_DEC_GFF if config["run"].get("emapper_kegganog") else []), 
        *(KEGGANOG if config["run"].get("emapper_kegganog") else []),
        QC_SUMMARY, 
        PANGENOME_PIRATE,
        PANGENOME_TREE, 
        PANGENOME_CORE_GENES_DISTANCE, 
        TYPING_MLST,
        TYPING_SPATYPER,
        ### Optional outputs depending on config cgmlst = true 
        *(CGMLST_PREP_EXTERNAL_DB if config["run"].get("cgmlst") else []), 
        *(CGMLST_ALLELE_CALL if config["run"].get("cgmlst") else []), 
        AMR_NCBI, 
        AMR_NCBI_SUMMARY,
        ### Optional outputs depending on config anvio = true 
        *(ANVIO_FASTA_REFORMATTED if config["run"].get("anvio") else []), 
        *(ANVIO_GENE_CALLS if config["run"].get("anvio") else []),
        *(ANVIO_ANNOTATION if config["run"].get("anvio") else []), 
        *(ANVIO_DB if config["run"].get("anvio") else []), 
        *(ANVIO_DONE_DB if config["run"].get("anvio") else []),
        *(ANVIO_EXTERNAL_GENOMES_FILE if config["run"].get("anvio") else []),
        *(ANVIO_STORAGE_DB if config["run"].get("anvio") else []), 
        *(ANVIO_PANGENOME_DB if config["run"].get("anvio") else []), 
        *(ANVIO_ANI if config["run"].get("anvio") else []), 
        ### Optional outputs depending on config snippy = true 
        *(VC_REFG_GBK if config["run"].get("snippy") else []),
        *(VC_SNIPPY if config["run"].get("snippy") else []),
        *(VC_COPY_VCF if config["run"].get("snippy") else []), 
        *(VC_COPY_VCF_GZ if config["run"].get("snippy") else []),
        *(VC_COPY_TAB if config["run"].get("snippy") else []), 
        *(VC_SNIPPY_CORE if config["run"].get("snippy") else []),
        *(VC_CORE_TREE if config["run"].get("snippy") else []),
        *(VC_VCF_HEATMAP if config["run"].get("snippy_vcf_heatmap") else []),
        *(VC_SNP_DISTANCE if config["run"].get("snippy") else []),
        ### Optional outputs depending on config gwas = true
        *(CSI if config["run"].get("gwas") else []), 
        *(VCF_MERGED if config["run"].get("gwas") else []),
        *(PHYLO_DIST if config["run"].get("gwas") else []),
        BAM,
        ### Optional outputs depending on config pirate_genes_of_interest = true
        *(PIRATE_GENE_TABLE if config["run"].get("filter_genes_of_interest") else []),
        *(SNP_GENES_OF_INTEREST if config["run"].get("filter_genes_of_interest") else []),
    message:
        "The genome analysis pipeline finished successfully!"
