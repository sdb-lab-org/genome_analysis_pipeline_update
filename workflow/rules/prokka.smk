###########################################
### Genome Annotation with Prokka        ##
###########################################

rule prokka:
    input:
        fixed_contigs       = OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.fasta"
    output:
        gff                 = OUTPUT_DIR + "{all_genomes}/annotation/prokka/prokka.gff",
        faa                 = OUTPUT_DIR + "{all_genomes}/annotation/prokka/prokka.faa",
        gbk                 = OUTPUT_DIR + "{all_genomes}/annotation/prokka/prokka.gbk",
        txt                 = OUTPUT_DIR + "{all_genomes}/annotation/prokka/prokka.txt",
        fna                 = OUTPUT_DIR + "{all_genomes}/annotation/prokka/prokka.fna",
        output_dir          = directory(OUTPUT_DIR + "{all_genomes}/annotation/prokka/")
    params:
        output_dir          = OUTPUT_DIR + "{all_genomes}/annotation/prokka/"
    message:
        "Running Prokka on fixed_contigs"
    conda: 
        "prokka_env"
    shell:
        "rm -rf {params.output_dir}; "
        "mkdir -p {params.output_dir}; "
        "prokka "
        "--outdir {params.output_dir} "
        "--prefix prokka "
        "--cpus {resources.cpus_per_task} "
        "{input.fixed_contigs} "
        "--force; "


rule copy_prokka_to_temp:
    input:
        faa                     = rules.prokka.output.faa,
        gff                     = rules.prokka.output.gff,
        txt                     = rules.prokka.output.txt
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
