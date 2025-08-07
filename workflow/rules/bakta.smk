###########################################
### Genome Annotation with Bakta        ###
###########################################

rule bakta:
    input:
        fixed_contigs = OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.fasta"
    output:
        gff           = OUTPUT_DIR + "{all_genomes}/annotation/bakta/bakta.gff3",
        faa           = OUTPUT_DIR + "{all_genomes}/annotation/bakta/bakta.faa",
        gbk           = OUTPUT_DIR + "{all_genomes}/annotation/bakta/bakta.gbff",
        txt           = OUTPUT_DIR + "{all_genomes}/annotation/bakta/bakta.txt"
    params:
        output_dir    = OUTPUT_DIR + "{all_genomes}/annotation/bakta/",
        bakta_db      = DATABASE_DIR + "bakta_db_full/"
    message:
        "Running Bakta on fixed_contigs"
    conda: 
        "bakta_env"
    shell:
        "rm -rf {params.output_dir}; "
        "mkdir -p {params.output_dir}; "
        "bakta  "
        "--output {params.output_dir} "
        "--prefix bakta "
        "--db {params.bakta_db} "
        "--threads {resources.cpus_per_task} "
        "{input.fixed_contigs} "
        "--verbose "
        "--force; "


rule copy_bakta_to_temp:
    input:
        faa                     = rules.bakta.output.faa,
        gff                     = rules.bakta.output.gff,
        txt                     = rules.bakta.output.txt
    output:
        temp_faa                = OUTPUT_DIR + "08_temp/temp_faa/{all_genomes}.faa",
        temp_gff                = OUTPUT_DIR + "08_temp/temp_gff/{all_genomes}.gff",
        temp_txt                = OUTPUT_DIR + "08_temp/temp_txt/{all_genomes}.txt"
    params:
        output_dir_faa          = OUTPUT_DIR + "08_temp/temp_faa/",
        output_dir_gff          = OUTPUT_DIR + "08_temp/temp_gff/",
        output_dir_txt          = OUTPUT_DIR + "08_temp/temp_txt/"
    message:
        "Copying annotaion files to temp folder"
    shell:
        "mkdir -p {params.output_dir_faa}; "
        "mkdir -p {params.output_dir_gff}; "
        "mkdir -p {params.output_dir_txt}; "
        "cp {input.faa} {output.temp_faa}; "
        "cp {input.gff} {output.temp_gff}; "
        "cp {input.txt} {output.temp_txt}; "
