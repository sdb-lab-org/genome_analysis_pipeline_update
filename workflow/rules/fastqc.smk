############################################
# Running FastQC on raw and filtered reads #
############################################
rule fastqc:
    input:
        read_1                  = WORKING_DIR + "raw_reads/{sample}_R1.fastq.gz",
        read_2                  = WORKING_DIR + "raw_reads/{sample}_R2.fastq.gz",
        trimmed_read_1          = rules.fastp.output.trimmed_1,
        trimmed_read_2          = rules.fastp.output.trimmed_2
    output:
        zip_file_1              = OUTPUT_DIR + "01_qc/fastqc/{sample}_R1_fastqc.zip",
        zip_file_2              = OUTPUT_DIR + "01_qc/fastqc/{sample}_R2_fastqc.zip",  
        trimmed_zip_file_1      = OUTPUT_DIR + "01_qc/fastqc/{sample}_trimmed_R1_fastqc.zip",
        trimmed_zip_file_2      = OUTPUT_DIR + "01_qc/fastqc/{sample}_trimmed_R2_fastqc.zip"        
    params:
        output_dir              = OUTPUT_DIR + "01_qc/fastqc/"     
    message:
        "Running FastQC on raw and trimmed files"
    conda: 
        "fastqc_env"
    shell:
        "echo {input.read_1}; "
        "echo {input.read_2}; "
        "echo {input.trimmed_read_1}; "
        "echo {input.trimmed_read_2}; "
        "mkdir -p {params.output_dir}; "
        "fastqc "
        "{input.read_1} {input.read_2} {input.trimmed_read_1} {input.trimmed_read_2} "
        "--outdir {params.output_dir} "
        "--threads {resources.cpus_per_task}"

