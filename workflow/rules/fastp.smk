############################################
## Running fastp to trimm reads           ##
############################################

rule fastp:
    input:
        read_1              = WORKING_DIR + "raw_reads/{sample}_R1.fastq.gz",
        read_2              = WORKING_DIR + "raw_reads/{sample}_R2.fastq.gz"
    output:   
        trimmed_1           = OUTPUT_DIR + "01_qc/fastp/{sample}_trimmed_R1.fastq.gz",
        trimmed_2           = OUTPUT_DIR + "01_qc/fastp/{sample}_trimmed_R2.fastq.gz",
        html                = OUTPUT_DIR + "01_qc/fastp/logs/{sample}_trimmed_fastp.html",
        json                = OUTPUT_DIR + "01_qc/fastp/logs/{sample}_trimmed_fastp.json"        
    message:
        "fastp trimming {wildcards.sample} reads"
    conda: 
        "fastp_env"
    params:
        cut_window_size     = config["fastp"]["cut_window_size"],
        cut_mean_quality    = config["fastp"]["cut_mean_quality"], 
        length_required     = config["fastp"]["length_required"],
        phread_quality      = config["fastp"]["phread_quality"],
        adapters            = config["refs"]["adapters"],
        dir                 = OUTPUT_DIR + "01_qc/fastp/"
    shell:
        "mkdir -p {params.dir}; "
        "fastp "
        "-i {input.read_1} "
        "-I {input.read_2} "
        "-o {output.trimmed_1} "
        "-O {output.trimmed_2} "
        "--thread {resources.cpus_per_task} "
        "--qualified_quality_phred {params.phread_quality} "
        "--cut_front "
        "--cut_tail "
        "--cut_window_size {params.cut_window_size} "
        "--cut_mean_quality {params.cut_mean_quality} "
        "--length_required {params.length_required} "
        "--html {output.html} "
        "--json {output.json} "
        "--detect_adapter_for_pe; "