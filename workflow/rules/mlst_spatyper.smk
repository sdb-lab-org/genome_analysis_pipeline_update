###########################################
### Typing                          ##
###########################################


rule run_mlst:
    input:
        fasta                   = expand(OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.fasta", all_genomes=ALL_GENOMES)
    output:
        mlst_tsv                = OUTPUT_DIR + "06_typing/mlst/mlst.tsv"
    params:
        fasta_dir               = OUTPUT_DIR + "08_temp/temp_fasta/",
        mlst_dir                = OUTPUT_DIR + "06_typing/mlst/"
    message:
        "Running MLST typing on all genome files"
    conda: 
        "mlst_env"
    shell:
        "mlst {params.fasta_dir}*.fasta > {output.mlst_tsv}"


rule spa_typing:
    input:
        contigs                 = OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.fasta"
    output:
        spatype                 = OUTPUT_DIR + "06_typing/spaTyper/{all_genomes}_spatype.txt"
    params:
        output_dir              = OUTPUT_DIR + "06_typing/spaTyper/"
    message:
        "Compute Spa Type Calculation"
    conda: 
        "spatyper_env"
    shell:
        "mkdir -p {params.output_dir}; "
        "spaTyper "
        "-f {input.contigs} "
        "--output {output.spatype}; "



