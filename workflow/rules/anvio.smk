
rule simplify_fasta_headers:
    input:
        fasta               = OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.fasta"
    output:
        simplified          = OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.simplified.fasta"
    conda:
        "base_env"  # or a specific environment if needed
    shell:
        "awk '/^>/ {{print $1}} !/^>/ {{print}}' {input.fasta} > {output.simplified}; "


rule parse_gff:
    input:
        gff                 = OUTPUT_DIR + "08_temp/temp_gff/{all_genomes}.gff"
    output:
        calls               = OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}_gene_calls.txt",
        annot               = OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}_gene_annot.txt"
    params:
        anvio_dir           = OUTPUT_DIR + "03_pangenome/anvio/db/"
    conda:
        "parser_env"
    shell:
        "mkdir -p {params.anvio_dir}; "
        "python scripts/gff_parser.py {input.gff} "
        "--gene-calls {output.calls} "
        "--annotation {output.annot}; "


rule create_contigs_db:
    input:
        fasta              = OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.simplified.fasta",
        gene_calls         = OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}_gene_calls.txt"
    output:
        db                 = OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}.db"
    params:
        project_name       = config["project_name"]
    conda:
        "anvio-8"
    shell:
        "anvi-gen-contigs-database "
        "-f {input.fasta} "
        "-o {output.db} "
        "--external-gene-calls {input.gene_calls} "
        "--project-name {params.project_name} "
        "--ignore-internal-stop-codons "
        "--num-threads {resources.cpus_per_task}; "


rule annotate_contigs_db:
    input:
        db                 = OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}.db",
        gene_annot         = OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}_gene_annot.txt"
    output:
        done               = OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}.annotation.done"
    conda:
        "anvio-8"
    shell:
        "anvi-import-functions "
        "-c {input.db} "
        "-i {input.gene_annot}; "
        "anvi-run-hmms "
        "-c {input.db}; "
        "touch {output.done}; "



rule create_external_genomes_file:
    input:
        dbs                         = expand(OUTPUT_DIR + "03_pangenome/anvio/db/{all_genomes}.db", all_genomes=ALL_GENOMES)
    output:
        external_genomes            = OUTPUT_DIR + "03_pangenome/anvio/external-genomes.txt"
    conda:
        "base"
    run:
        os.makedirs(os.path.dirname(output.external_genomes), exist_ok=True)
        with open(output.external_genomes, "w") as f:
            f.write("name\tcontigs_db_path\n")
            for path in input.dbs:
                raw_sample = os.path.basename(path).replace(".db", "")
                # Replace disallowed characters with underscores
                sanitized_sample = re.sub(r"[^A-Za-z0-9_]", "_", raw_sample)
                f.write(f"{sanitized_sample}\t{os.path.abspath(path)}\n")


rule create_genomes_storage:
    input:
        external_genomes        = OUTPUT_DIR + "03_pangenome/anvio/external-genomes.txt"
    output:
        storage                 = OUTPUT_DIR + "03_pangenome/anvio/anvio_storage-GENOMES.db"
    conda:
        "anvio-8"
    shell:
        "anvi-gen-genomes-storage "
        "-e {input.external_genomes} "
        "-o {output.storage} "
        "--gene-caller Prodigal; "


rule run_anvio_pangenome:
    input:
        storage         = OUTPUT_DIR + "03_pangenome/anvio/anvio_storage-GENOMES.db"
    output:
        pan_db_txt      = OUTPUT_DIR + "03_pangenome/anvio/anvio_Pangenome-PAN.txt"
    params:
        output_dir      = OUTPUT_DIR + "03_pangenome/anvio/"
    conda:
        "anvio-8"
    shell:
        "anvi-pan-genome "
        "-g {input.storage} "
        "--project-name anvio_Pangenome "
        "--output-dir {params.output_dir} "
        "--num-threads {resources.cpus_per_task} "
        "--minbit 0.5 "
        "--mcl-inflation 10 "
        "--use-ncbi-blast; "
        "touch {output.pan_db_txt}; "


rule run_ani_similarity:
    input:
        external_genomes    = OUTPUT_DIR + "03_pangenome/anvio/external-genomes.txt",
        pan_db_txt      = OUTPUT_DIR + "03_pangenome/anvio/anvio_Pangenome-PAN.txt"
    output:
        ani_txt             = OUTPUT_DIR + "03_pangenome/anvio/ANI/genome_similarity.txt"
    params:
        ani_dir             = OUTPUT_DIR + "03_pangenome/anvio/ANI/",
        anvio_dir           = OUTPUT_DIR + "03_pangenome/anvio/"
    conda:
        "anvio-8"
    shell:
        "mkdir -p {params.ani_dir}; "
        "anvi-compute-genome-similarity "
        "--external-genomes {input.external_genomes} "
        "--program pyANI "
        "--output-dir {params.ani_dir} "
        "--force-overwrite "
        "--num-threads {resources.cpus_per_task} "
        "--pan-db {params.anvio_dir}anvio_Pangenome-PAN.db; "
        "touch {output.ani_txt}; "
