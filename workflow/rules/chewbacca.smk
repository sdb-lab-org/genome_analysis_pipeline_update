rule chewbbaca_prep_schema:
    input:
        schema_dir          = CGMLST_SCHEMA_DIR + "alleles"
    output:
        prepared_schema     = CGMLST_SCHEMA_DIR + "PreparedSchema/schema.json"
    params:
        output_dir          = CGMLST_SCHEMA_DIR + "PreparedSchema"
    conda:
        "chewie_env"
    message:
        "Preparing chewBBACA schema allele FASTAs"
    shell:
        "rm -rf {params.output_dir}; " 
        "chewBBACA.py PrepExternalSchema "
        "-i {input.schema_dir} "
        "-o {params.output_dir}; "
        "touch {output.prepared_schema}; "




rule chewbbaca_allele_call:
    input:
        data            = expand(OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.fasta", all_genomes=ALL_GENOMES),
        edited_schema   = CGMLST_SCHEMA_DIR + "PreparedSchema"
    output:
        allele_profile  = OUTPUT_DIR + "06_typing/chewbbaca/allele_call_out/results_alleles.tsv",
        allele_stats    = OUTPUT_DIR + "06_typing/chewbbaca/allele_call_out/results_statistics.tsv",
        contig_info     = OUTPUT_DIR + "06_typing/chewbbaca/allele_call_out/results_contigsInfo.tsv",
        distances       = OUTPUT_DIR + "06_typing/chewbbaca/distances.tab",
        alleles_call    = OUTPUT_DIR + "06_typing/chewbbaca/results_allele_call.json",
        cgmlst_profile   = OUTPUT_DIR + "06_typing/chewbbaca/extract_cgmlst/cgMLST.tsv"
    params:
        data            = OUTPUT_DIR + "08_temp/temp_fasta/",
        outdir          = OUTPUT_DIR + "06_typing/chewbbaca",
        extract_out     = OUTPUT_DIR + "06_typing/chewbbaca/extract_cgmlst",
        tmp_call_dir    = OUTPUT_DIR + "06_typing/chewbbaca/tmp_allele_call"
    conda:
        "chewie_env"
    shell:
        """
        rm -rf {params.tmp_call_dir}
        rm -rf {params.outdir}/allele_call_evaluation

        chewBBACA.py AlleleCall -i {params.data} -g {input.edited_schema} -o {params.tmp_call_dir} --cpu {resources.cpus_per_task}
        chewBBACA.py AlleleCallEvaluator --input-files {params.tmp_call_dir} -g {input.edited_schema} -o {params.outdir}/allele_call_evaluation
        chewBBACA.py ExtractCgMLST -i {params.tmp_call_dir}/results_alleles.tsv -o {params.extract_out}
        cgmlst-dists {params.tmp_call_dir}/results_alleles.tsv > {output.distances}

        cp {params.tmp_call_dir}/results_alleles.tsv {output.allele_profile}
        cp {params.tmp_call_dir}/results_statistics.tsv {output.allele_stats}
        cp {params.tmp_call_dir}/results_contigsInfo.tsv {output.contig_info}

        touch {output.alleles_call}
        """
