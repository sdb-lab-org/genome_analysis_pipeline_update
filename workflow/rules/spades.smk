##################################################
## Genome Assembly with SPAdes                  ##
##################################################

rule spades:
    input:
        trimmed_1           = rules.fastp.output.trimmed_1,
        trimmed_2           = rules.fastp.output.trimmed_2
    output:
        contigs_spades      = OUTPUT_DIR + "{sample}/assembly/spades/contigs.fasta"
    message:
        "Run Spades to assemble trimmed reads"
    conda: 
        "spades_env"
    params:
        output_dir          = OUTPUT_DIR + "{sample}/assembly/spades/"
    shell:
        "mkdir -p {params.output_dir}; "
        "spades.py "
        "-1 {input.trimmed_1} "
        "-2 {input.trimmed_2} "
        "--careful "
        "--threads {resources.cpus_per_task} " 
        "-o {params.output_dir}; "



rule fix_contigs_spades:
    input:
        contigs             = rules.spades.output.contigs_spades
    output:
        fixed_contigs       = OUTPUT_DIR + "{sample}/assembly/spades/contigs_fixed.fasta"
    message:
        "Fix lenght of contig names and remove contigs < 200bp"
    conda: 
        "seqkit_env"
    shell:
        """
        seqkit seq -m 200 {input.contigs} | awk '/^>/ {{print substr($0,1,20); next}} 1' > {output.fixed_contigs}
        """


rule copy_spades_to_temp:
    input:
        fasta                   = rules.fix_contigs_spades.output.fixed_contigs
    output:
        temp_fasta              = OUTPUT_DIR + "08_temp/temp_fasta/{sample}.fasta"
    params:
        output_dir_fasta        = OUTPUT_DIR + "08_temp/temp_fasta/"
    message:
        "Copying annotaion files to temp folder"
    shell:
        "mkdir -p {params.output_dir_fasta}; "
        "cp {input.fasta} {output.temp_fasta}; "



