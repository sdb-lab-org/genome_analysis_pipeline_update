
rule align_reads_to_reference:
    input:
        ref             = lambda wildcards: f"{OUTPUT_DIR}08_temp/temp_fasta/{REF_SAMPLE}.fasta",
        reads1          = rules.fastp.output.trimmed_1,
        reads2          = rules.fastp.output.trimmed_2
    output:
        bam             = OUTPUT_DIR + "{sample}/mapping/{sample}.sorted.bam",
        bai             = OUTPUT_DIR + "{sample}/mapping/{sample}.sorted.bam.bai"
    params:
        output_dir      = OUTPUT_DIR + "{sample}/mapping/"
    conda:
        "mapping_env"  
    shell:
        """
        bwa index {input.ref}
        mkdir -p {params.output_dir}

        bwa mem -t {resources.cpus_per_task} {input.ref} {input.reads1} {input.reads2} |
        samtools view -bS - |
        samtools sort -@ {resources.cpus_per_task} -o {output.bam}

        samtools index {output.bam}
        """
