###################################################
## Variant calling and visualization             ##
###################################################

rule snippy:
    input:
        fixed_contigs       = OUTPUT_DIR + "08_temp/temp_fasta/{all_genomes}.fasta",
        refg                = rules.refg_gbk.output.gbk
    output:
        snippy_output       = OUTPUT_DIR + "04_variant_calling/snippy/{all_genomes}/snps.vcf",
        snippy_output_gz    = OUTPUT_DIR + "04_variant_calling/snippy/{all_genomes}/snps.vcf.gz",
        snippy_tab          = OUTPUT_DIR + "04_variant_calling/snippy/{all_genomes}/snps.tab",
    params:
        output_dir          = OUTPUT_DIR + "04_variant_calling/snippy/{all_genomes}/"
    conda: "snippy_env"
    shell:
        "mkdir -p {params.output_dir}; "
        "snippy "
        "--report "
        "--cpus {resources.cpus_per_task} "
        "--outdir {params.output_dir} "
        "--ref {input.refg} "
        "--ctgs {input.fixed_contigs} "
        "--force; "

rule copy_snippy_to_temp:
    input:
        vcf                     = rules.snippy.output.snippy_output,
        vcf_gz                  = rules.snippy.output.snippy_output_gz,
        tab                     = rules.snippy.output.snippy_tab
    output:
        temp_vcf                = OUTPUT_DIR + "08_temp/temp_vcf/{all_genomes}.snps.vcf",
        temp_vcf_gz             = OUTPUT_DIR + "08_temp/temp_vcf_gz/{all_genomes}.snps.vcf.gz",
        temp_tab                = OUTPUT_DIR + "08_temp/temp_tab/{all_genomes}.snps.tab"
    params:
        output_dir_vcf          = OUTPUT_DIR + "08_temp/temp_vcf/",
        output_dir_vcf_gz       = OUTPUT_DIR + "08_temp/temp_vcf_gz/",
        output_dir_tab          = OUTPUT_DIR + "08_temp/temp_tab/"
    message:
        "Copying snippy vcf files to temp folder"
    shell:
        "mkdir -p {params.output_dir_vcf}; "
        "mkdir -p {params.output_dir_vcf_gz}; "
        "mkdir -p {params.output_dir_tab}; "
        "echo {input.vcf}; "
        "cp {input.vcf} {output.temp_vcf}; "
        "cp {input.vcf_gz} {output.temp_vcf_gz}; "
        "cp {input.tab} {output.temp_tab}; "




rule snippy_core:
    input:
        refg                       = expand(rules.refg_gbk.output.gbk, all_genomes=ALL_GENOMES)
    output:
        aln                        = OUTPUT_DIR + "04_variant_calling/snippy-core/core.aln",
        vcf                        = OUTPUT_DIR + "04_variant_calling/snippy-core/core.vcf",
        full_aln                   = OUTPUT_DIR + "04_variant_calling/snippy-core/core.full.aln"
    params:
        prefix                      = OUTPUT_DIR + "04_variant_calling/snippy-core/core",
        snippy_dirs                 = lambda wildcards: " ".join([
            OUTPUT_DIR + f"04_variant_calling/snippy/{all_genomes}"
            for all_genomes in ALL_GENOMES if all_genomes != REF_SAMPLE
        ]),
        output_dir                  = OUTPUT_DIR + "04_variant_calling/snippy-core/"
    message:
        "Aligning core and whole genomes into a multi fasta file"
    conda: 
        "snippy_env"
    shell:
        "echo 'Reference genome: {input.refg}'; "
        "echo 'Including snippy directories: {params.snippy_dirs}'; "
        "mkdir -p {params.output_dir}; "
        "snippy-core "
        "--ref {input.refg} "
        "--prefix {params.prefix} "
        "{params.snippy_dirs}; "

        
rule tree:
    input:
        aln                     = rules.snippy_core.output.full_aln
    output:
        iqtree_log              = OUTPUT_DIR + "04_variant_calling/snippy-core/iqtree.log"
    params:
        dir                     = OUTPUT_DIR + "04_variant_calling/snippy-core/"
    message:
        "Builing phylogeny tree of whole genomes using IQ-Tree"
    conda: 
        "iqtree_env"
    shell:
        "iqtree "
        "-s {input.aln} "
        "-bb 1000 "
        "> {output.iqtree_log}; "

rule vcf_viewer:
    input:
        vcf                     = rules.snippy_core.output.vcf
    output:
        heatmap_figure          = OUTPUT_DIR + "04_variant_calling/vcf_viewer/heatmap_output.html"
    params:
        output_dir              = OUTPUT_DIR + "04_variant_calling/vcf_viewer/"
    message:
        "Generating SNP heatmap"
    shell:
        "mkdir -p {params.output_dir}; "
        "Rscript scripts/vcf2heatmap.R {input.vcf} > {output.heatmap_figure}"



rule snp_dists:
    input:
        aln                     = rules.snippy_core.output.full_aln
    output:
        dist_matrix             = OUTPUT_DIR + "04_variant_calling/snippy-core/snp_distance_matrix.tsv"
    message:
        "Calculating SNP distance matrix from core genome alignment (snippy-core)"
    conda:
        "snp-dists_env"
    shell:
        "snp-dists {input.aln} > {output.dist_matrix}"


    
