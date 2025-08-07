############################################
## Running quast for qc of assemblies     ##
############################################

rule quast:
    input:
        fixed_contigs           = OUTPUT_DIR + "08_temp/temp_fasta/{sample}.fasta"
    output:
        report                  = OUTPUT_DIR + "01_qc/quast/{sample}/report.txt"
    params:
        output_dir              = OUTPUT_DIR + "01_qc/quast/{sample}/"
    message:
        "Running Quast on fixed_contigs file"
    conda: 
        "quast_env"
    shell:
        "mkdir -p {params.output_dir}; "
        "quast.py "
        "-o {params.output_dir} "
        "--threads {resources.cpus_per_task} "
        "{input.fixed_contigs}"
