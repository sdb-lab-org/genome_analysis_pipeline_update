############################################
## MultiQC report                         ##
############################################

rule multiqc:
    input:
        zip_file_1              = expand(rules.fastqc.output.zip_file_1, sample=SAMPLES),
        zip_file_2              = expand(rules.fastqc.output.zip_file_2, sample=SAMPLES),
        trimmed_zip_file_1      = expand(rules.fastqc.output.trimmed_zip_file_1, sample=SAMPLES),
        trimmed_zip_file_2      = expand(rules.fastqc.output.trimmed_zip_file_2, sample=SAMPLES)  
    output:
        output                  = OUTPUT_DIR + "01_qc/multiqc/multiqc_report.html"
    params:
        dir_fastqc              = OUTPUT_DIR + "01_qc/fastqc/",
        outdir                  = OUTPUT_DIR + "01_qc/multiqc/"
    message: 
        "Summarising reports with multiqc"
    conda: 
        "multiqc_env"
    shell:
        "mkdir -p {params.outdir}; "
        "multiqc "
        "--force "
        "--outdir {params.outdir} "
        "{params.dir_fastqc}; "
