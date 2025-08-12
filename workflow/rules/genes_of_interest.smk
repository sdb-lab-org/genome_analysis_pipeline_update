rule pirate_genes_of_interest:
    input:
        pirate_table                = OUTPUT_DIR + "03_pangenome/pirate/PIRATE.gene_families.ordered.tsv"
    output:
        presence_absence_table      = OUTPUT_DIR + "03_pangenome/genes_of_interest_presence_absence.tsv"
    message:
        "Compare genome content using PIRATE and filter for genes of interest"
    run:
        # Load PIRATE gene families table
        df = pd.read_csv(input.pirate_table, sep="\t")

        # Load genes of interest from config
        genes_of_interest = config["genes_of_interest"]

        # The genome columns (everything after the PIRATE metadata columns)
        meta_cols = [
            "allele_name", "gene_family", "consensus_gene_name", "consensus_product",
            "threshold", "alleles_at_maximum_threshold", "number_genomes",
            "average_dose", "min_dose", "max_dose",
            "genomes_containing_fissions", "genomes_containing_duplications",
            "number_fission_loci", "number_duplicated_loci", "no_loci",
            "products", "gene_names", "min_length(bp)", "max_length(bp)",
            "average_length(bp)", "cluster", "cluster_order"
        ]
        genome_cols = [c for c in df.columns if c not in meta_cols]

        # Filter to genes of interest
        filtered_df = df[df["consensus_gene_name"].isin(genes_of_interest)]

        # Build presence/absence table
        presence_absence = []
        for gene in genes_of_interest:
            row = {"genes_of_interest": gene}
            rows_for_gene = filtered_df[filtered_df["consensus_gene_name"] == gene]
            for genome in genome_cols:
                # presence if ANY non-null/non-empty value in this genome column
                present = (
                    rows_for_gene[genome].notna() &
                    (rows_for_gene[genome].astype(str).str.strip() != "")
                ).any()
                row[genome] = 1 if present else 0
            presence_absence.append(row)

        pa_df = pd.DataFrame(presence_absence, columns=["genes_of_interest"] + genome_cols)

        # Write to file
        pa_df.to_csv(output.presence_absence_table, sep="\t", index=False)


rule snp_genes_of_interest:
    input:
        snptab                      = expand(OUTPUT_DIR + "08_temp/temp_tab/{all_genomes}.snps.tab", all_genomes=ALL_GENOMES)
    output:
        matrix                      = OUTPUT_DIR + "04_variant_calling/genes_of_interest_snp_matrix.tsv"
    params:
        output_dir                  = OUTPUT_DIR + "04_variant_calling/",
        genes                       = lambda wildcards: " ".join(config["genes_of_interest"])
    conda:
        "r_env"
    message:
        "Building SNP presence/absence matrix (with EFFECT) for selected genes: {params.genes}"
    shell:
        """
        mkdir -p {params.output_dir}
        Rscript workflow/scripts/snptab2snp_matrix.R {input.snptab} "{params.genes}" {output.matrix}
        """
