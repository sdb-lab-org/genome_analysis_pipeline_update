###################################################
### AMR/virulence genes screening with ABRICATE  ##
###################################################

rule abricate:
    input:
        fasta                   = OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.fasta"
    output:
        ncbi                    = OUTPUT_DIR + "07_amr/{all_genomes}.ncbi.tsv",
        resfinder               = OUTPUT_DIR + "07_amr/{all_genomes}.resfinder.tsv",
        plasmidfinder           = OUTPUT_DIR + "07_amr/{all_genomes}.plasmidfinder.tsv",
        card                    = OUTPUT_DIR + "07_amr/{all_genomes}.card.tsv",
        vfdb                    = OUTPUT_DIR + "07_amr/{all_genomes}.vfdb.tsv"
    params:
        output_dir              = OUTPUT_DIR + "07_amr/"
    message:
        "Screen for amr genes"
    conda: 
        "abricate_env"
    shell:
        "mkdir -p {params.output_dir}; "
        "abricate --db ncbi {input.fasta} --nopath > {output.ncbi}; "
        "abricate --db resfinder {input.fasta} --nopath > {output.resfinder}; "
        "abricate --db plasmidfinder {input.fasta} --nopath > {output.plasmidfinder}; "
        "abricate --db card {input.fasta} --nopath > {output.card}; "
        "abricate --db vfdb {input.fasta} --nopath > {output.vfdb}; "

rule abricate_table:
    input:
        ncbi                    = expand(OUTPUT_DIR + "07_amr/{all_genomes}.ncbi.tsv", all_genomes=ALL_GENOMES),
        resfinder               = expand(OUTPUT_DIR + "07_amr/{all_genomes}.resfinder.tsv", all_genomes=ALL_GENOMES),
        plasmidfinder           = expand(OUTPUT_DIR + "07_amr/{all_genomes}.plasmidfinder.tsv", all_genomes=ALL_GENOMES),
        card                    = expand(OUTPUT_DIR + "07_amr/{all_genomes}.card.tsv", all_genomes=ALL_GENOMES),
        vfdb                    = expand(OUTPUT_DIR + "07_amr/{all_genomes}.vfdb.tsv", all_genomes=ALL_GENOMES)
    output:
        summary_ncbi            = OUTPUT_DIR + "07_amr/amr_ncbi_summary.tab",
        summary_resfinder       = OUTPUT_DIR + "07_amr/amr_resfinder_summary.tab",
        summary_plasmidfinder   = OUTPUT_DIR + "07_amr/amr_plasmidfinder_summary.tab",
        summary_card            = OUTPUT_DIR + "07_amr/amr_card_summary.tab",
        summary_vfdb            = OUTPUT_DIR + "07_amr/amr_vfdb_summary.tab"
    params:
        output_dir              = OUTPUT_DIR + "07_amr/"
    message:
        "Screen for amr genes"
    conda: 
        "abricate_env"
    shell:
        "abricate --summary {params.output_dir}*.ncbi.tsv > {output.summary_ncbi}; "
        "abricate --summary {params.output_dir}*.resfinder.tsv > {output.summary_resfinder}; "
        "abricate --summary {params.output_dir}*.plasmidfinder.tsv > {output.summary_plasmidfinder}; "
        "abricate --summary {params.output_dir}*.card.tsv > {output.summary_card}; "
        "abricate --summary {params.output_dir}*.vfdb.tsv > {output.summary_vfdb}; "




