rule copy_genome:
    input:
        fasta           = WORKING_DIR + "genomes/{genome}.fasta"
    output:
        contigs         = OUTPUT_DIR + "08_temp/temp_fasta/{genome}.fasta"
    params:
        outdir          = WORKING_DIR + "genomes/"
    shell:
        "cp {input.fasta} {output.contigs}; "
