###################################################
### Compare genome content (PIRATE)              ##
###################################################

rule pirate:
    input:
        gff_files           = expand(OUTPUT_DIR + "08_temp/temp_gff/{all_genomes}.gff", all_genomes=ALL_GENOMES)
    output:
        alignment_fasta     = OUTPUT_DIR + "03_pangenome/pirate/core_alignment.fasta",
        gene_families       = OUTPUT_DIR + "03_pangenome/pirate/PIRATE.gene_families.ordered.tsv" 
    params:
        gff_dir             = OUTPUT_DIR + "08_temp/temp_gff/",
        output_dir          = OUTPUT_DIR + "03_pangenome/pirate/"
    message:
        "Compaire genome content using pirate"
    conda: 
        "pirate_env"
    shell:
        "mkdir -p {params.output_dir}; "
        "PIRATE "
        "-i {params.gff_dir} "
        "-s 50,70,90,95,98 "
        "-o {params.output_dir} "
        "-t {resources.cpus_per_task} "
        "-a "
        "-r; "



rule pirate_prevlocustag:
    input:
        gene_families = OUTPUT_DIR + "03_pangenome/pirate/PIRATE.gene_families.ordered.tsv"
    output:
        gene_families_prev = OUTPUT_DIR + "03_pangenome/pirate/PIRATE.gene_families.ordered_prevlocustag.tsv"
    params:
        modified_gffs = OUTPUT_DIR + "03_pangenome/pirate/modified_gffs/"
    conda:
        "pirate_env"   
    shell:
        "perl workflow/scripts/subsample_outputs.pl "
        "-i {input.gene_families} "
        "-g {params.modified_gffs} "
        "-o {output.gene_families_prev} "
        "--field prev_locus;"

rule tree_pirate:
    input:
        aln                     = rules.pirate.output.alignment_fasta
    output:
        iqtree_log              = OUTPUT_DIR + "03_pangenome/pirate/iqtree.log",
        treefile                = OUTPUT_DIR + "03_pangenome/pirate/core_alignment.fasta.treefile"
    params:
        dir                     = OUTPUT_DIR + "03_pangenome/pirate/"
    message:
        "Builing phylogeny tree of whole genomes using IQ-Tree"
    conda: 
        "iqtree_env"
    shell:
        "iqtree "
        "-s {input.aln} "
        "-bb 1000 "
        "> {output.iqtree_log}; "



rule core_gene_snp_dists:
    input:
        aln                     = rules.pirate.output.alignment_fasta
    output:
        dist_matrix             = OUTPUT_DIR + "03_pangenome/pirate/core_gene_snp_distances.tsv"
    message:
        "🔍 Calculating SNP distances between genomes based on core gene alignments (from PIRATE output)"
    conda:
        "snp-dists_env"
    shell:
        "snp-dists {input.aln} > {output.dist_matrix}"
