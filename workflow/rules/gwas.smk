rule bcf_index:
    input:
        vcf                  = OUTPUT_DIR + "08_temp/temp_vcf_gz/{gwas_genomes}.snps.vcf.gz"
    output:
        csi                  = OUTPUT_DIR + "08_temp/temp_vcf_gz/{gwas_genomes}.snps.vcf.gz.csi"
    params:
        dir                  = OUTPUT_DIR + "05_gwas/"
    message:
        "Genome wide association study"
    conda:
        "pyseer_env"
    shell:
        "mkdir -p {params.dir}; "
        "bcftools index "
        "-f "
        "-o {output.csi} "
        "{input.vcf}; "

rule vcf_merge:
    input:
        vcf                  = expand(OUTPUT_DIR + "08_temp/temp_vcf_gz/{gwas_genomes}.snps.vcf.gz", gwas_genomes=GWAS_GENOMES) 
    output:
        merged_vcf           = OUTPUT_DIR + "05_gwas/merged.vcf.gz"
    params:
        dir                  = OUTPUT_DIR + "05_gwas/"
    message:
        "Genome wide association study"
    conda:
        "pyseer_env"
    shell:
        "mkdir -p {params.dir}; "
        "bcftools merge {input.vcf} "
        "-Oz "
        "-o {output.merged_vcf}; "


rule phylogeny_distance:
    input:
        tree            = OUTPUT_DIR + "03_pangenome/pirate/core_alignment.fasta.treefile" 
    output:
        tsv             = OUTPUT_DIR + "05_gwas/phylogeny_dists.tsv"
    conda:
        "pyseer_env"  
    shell:
        "python scripts/phylogeny_distance.py {input.tree} > {output.tsv}"

