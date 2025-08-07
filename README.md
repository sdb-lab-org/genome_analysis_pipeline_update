# Genome Analysis Pipeline Update

This Snakemake-based pipeline automates the analysis of paired-end sequencing data for microbial genomes. It performs genome assembly, functional annotation, variant calling, pangenome analysis, and typing (e.g., spa typing, MLST, cgMLST), along with antimicrobial resistance (AMR) and virulence gene screening.

## Overview
The Snakemake pipeline performs the following steps (by default/optional):

**Quality Control and Trimming (default)**
- `fastp`: Performs quality filtering and adapter trimming on raw sequencing reads using [fastp](https://github.com/OpenGene/fastp).
- `fastqc`: Runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on both raw and trimmed reads to assess read quality.
- `multiqc`: Aggregates all FastQC reports into a single interactive summary using [MultiQC](https://multiqc.info/).

**Genome Assembly (default)**
- `spades`: Assembles trimmed reads into contigs using [SPAdes](http://cab.spbu.ru/software/spades/) and contigs <= 200 bp are filtered out (default).
- `unicycler`: Alternatively, runs [Unicycler](https://github.com/rrwick/Unicycler) if specified in config.yaml (assembler: unicycler).
- `quast`: Evaluates assembly quality using [QUAST](http://quast.sourceforge.net/quast).

**Genome Annotation and Functional Prediction**
- `prokka`:  Annotates assembled contigs with [Prokka](https://github.com/tseemann/prokka), identifying genes, rRNAs, tRNAs, and other features (default).
- `bakta`: Annotates assembled contigs with [Bakta](https://github.com/oschwengers/bakta) if specified in config.yaml (annotator: bakta).
- `emapper_kegganog`: Maps predicted proteins to orthologous groups and functional categories using [eggNOG-mapper](http://eggnog-mapper.embl.de/) and integrates KEGG pathway data with eggNOG annotations to provide insights into metabolic and functional pathways using [KEGGaNOG](https://github.com/iliapopov17/KEGGaNOG) (optional).
   	  
**Pangenome Analysis**
- `pirate`: Runs [PIRATE](https://github.com/SionBayliss/PIRATE) for pangenome analysis, identifying core and accessory genes across multiple genomes (default).
- `iqtree`: Generates a tree out of the core SNP alignment using [IQ-TREE](http://www.iqtree.org) (default)
- `snp-dists`: Calculates pairwise nucleotide differences from the core gene-by-gene alignment using [snp-dists](https://github.com/tseemann/snp-dists) (default).
- `anvio`: Runs [Anvi’o](https://anvio.org) pangenomic analysis and creates databases for a ringplot (optional) (-->see some solutions). 
  
**Variant Calling and Visualization**
- `snippy`: Detects single nucleotide polymorphisms (SNPs) relative to the reference genome using [Snippy](https://github.com/tseemann/snippy) (default).
- `snippy-core`: Combines the SNPs from multiple samples to create a core SNP alignment for phylogenetic analysis, generates a tree out of the core SNP alignment using [IQ-TREE](http://www.iqtree.org) (default).
- `snp-dists`: Calculates pairwise SNP differences from the core SNP alignment using [snp-dists](https://github.com/tseemann/snp-dists) (default).
- `vcf_viewer`: Generates a heatmap to visualize variations across strains (optional) --> difficult to execute if there are too many snps, but you could filter them first?.
  
**AMR/virulence genes screening**
- `abricate`: Screens genomes for antimicrobial resistance and virulence genes using [ABRicate](https://github.com/tseemann/abricate) blasting against the virulence factor database [VFDB](http://www.mgc.ac.cn/VFs/), [NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047), [CARD](https://card.mcmaster.ca), [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder) and [PlasmidFinder](https://cge.cbs.dtu.dk/services/PlasmidFinder). Results from each database are summarized individually to be visualized as a heatmap (default).
   
**Typing**
- `spa_typing`: Uses [spaTyper](https://github.com/medvir/spaTyper) to determine spa types from the assembled contigs for characterizing Staphylococcus aureus strains (default).
- `mlst`: Uses [mlst]([https://github.com/tseemann/mlst)]) to determine the MLST types from the assembled contigs (default).
- `chewbacca`: Uses [chewBBACA](https://github.com/B-UMMI/chewBBACA) to determine the cgMLST of the assembled contigs--> see soem solutions (optional).

**Filtering for genes of interest (optional) -->genes must be annotated in your reference genome (e.g., gene="tpiA")**
- `pirate_genes_of_interest`: Extracts the genes specified in the config file (genes_of_interest) from the PIRATE output and generates a presence/absence table across all sample.
- `snp_genes_of_interest`: Searches for all SNPs within the genes specified in the config file and creates a presence/absence matrix of these SNP positions across all samples.
Note: Only SNPs located in annotated genes are considered.

## Directory Structure
```
├── config
│   └── config.yaml               # Config file for parameters & paths  
├── inputs
│   ├── adapters
│   │   └── adapters.fa           # Adapter sequences for trimming
│   ├── genomes			  # External genomes (See some solutions section)    
│   ├── raw_reads                 # Input FASTQ files
│   │   ├── C22_R1.fastq.gz
│   │   ├── C22_R2.fastq.gz
│   │   ├── C24_R1.fastq.gz
│   │   └── C24_R2.fastq.gz
│   └── reference
│       └── refg.fasta            # Reference genome (optional)   
├── output                        # All results
└── workflow                      # Snakemake rules, envs, scripts
    ├── envs
    ├── log
    ├── profiles
    │   ├── default
    │   └── ga_pipeline
    ├── rules
    └── scripts
```


## Installation

1. Clone this git repository
   ```bash
   git clone https://github.com/doina1234/genome_analysis_pipeline_update.git
   ```
   
2. If the conda environments are not installed on your computer, install using `conda env create -f <environment>.yml` command. The environment files are in `workflow/envs` directory.
   ```bash
   conda env create -f workflow/envs/<environment>.yml
   ```
   
3. Eggnog database: Follow setup instructions eggnog-mapper documentation (https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12). I recommend creating a databases folder and adding the eggnog-mapper-data folder in there. After a successful download, add the path to the databases in the `config.yaml file`, something like `path/databases/eggnog-mapper-data`.

4. Bakta database: If you want to use the bakta annotation tool (config.yaml: annotator: bakta), follow the bakta database setup instructions (https://github.com/oschwengers/bakta?tab=readme-ov-file#database) and add the path to the databases in the `config.yaml file`, something like `path/databases/bakta_db_full`.
   
5. Anvi’o: It is recommended to test the installation of all required tools, but especially verifying that Anvi’o is correctly installed and functional by running its built-in test suites:
   ```bash
   conda activate anvio-8
   anvi-self-test --suite mini
   anvi-self-test --suite pangenomics
   conda deactivate
   ```


## Setup

1. Prepare the `config.yaml` file in the `config/` directory.
    ```
   ########################################
   ## Configurations                     ##
   ########################################
   
   # Define samples:
   samples:
     - {sample}
   genomes:
     - {genome}
   refg: 				# if you provide a reference genome (optional, put the .fasta file into the ../inputs/genomes/ folder 					# and write here the sample_name (eg - C22)
     - {reference}
   
   # Define additional rules to execute:
   run:
     emapper_kegganog: false       	# if you want to run eggnog-mapper and KEGG annotation
     snippy: true                  	# if you want to run snippy
     snippy_vcf_heatmap: false     	# option: true if you want to generate a heatmap from snippy-core vcf file
     gwas: false                   	# if you want to run GWAS analysis 
     anvio: false                  	# if you want to run anvio 
     cgmlst: false                 	# if you want to run cgmlst analysis
     filter_genes_of_interest: false  	# if you want to filter genes of interest

   # Define the gene of interest:
   genes_of_interest:
     - {gene}
   
   # Define the directories:
   output_dir: output/
   working_dir: inputs/  
   database_dir: /path_to/database/
   project_dir: /path_to/projects/{project_name}/
   project_name: {project_name}
   cgmlst_schema: /path_to/database/{cgmlst_allele_folder}
   
   # Define annotator (default prokka, other options: bakta)
   annotator: prokka
   
   # Define assembler (default spades, other options: unicycler)
   assembler: spades
   
   # Define parameters for specific tools: 
   refs:
     adapters: inputs/adapters/adapters.fa
   
   fastp:
     phread_quality: 20      # Phred+33 score
     cut_window_size: 4      # 4bp sliding window
     cut_mean_quality: 15    # Trim when window average < 20
     length_required: 36     # Discard reads shorter than 50 bp
   ```    ```
   ########################################
   ## Configurations                     ##
   ########################################
   
   # Define samples:
   samples:
     - {sample}
   genomes:
     - {genome}
   refg: 
     - {reference}
   
   # Define additional rules to execute:
   run:
     emapper_kegganog: false       # if you want to run eggnog-mapper and KEGG annotation
     snippy: true                  # if you want to run snippy
     snippy_vcf_heatmap: false     # option: true if you want to generate a heatmap from snippy-core vcf file
     gwas: false                   # if you want to run GWAS analysis 
     anvio: false                  # if you want to run anvio 
     cgmlst: false                 # if you want to run cgmlst analysis
     filter_genes_of_interest: true  # if you want to filter genes of interest

   # Define the gene of interest:
   genes_of_interest:
     - {gene}
   
   # Define the directories:
   output_dir: output/
   working_dir: inputs/  
   database_dir: /path_to/database/
   project_dir: /path_to/projects/{project_name}/
   project_name: {project_name}
   cgmlst_schema: /path_to/database/{cgmlst_allele_folder}
   
   # Define annotator (default prokka, other options: bakta)
   annotator: prokka
   
   # Define assembler (default spades, other options: unicycler)
   assembler: spades
   
   # Define parameters for specific tools: 
   refs:
     adapters: inputs/adapters/adapters.fa
   
   fastp:
     phread_quality: 20      # Phred+33 score
     cut_window_size: 4      # 4bp sliding window
     cut_mean_quality: 15    # Trim when window average < 20
     length_required: 36     # Discard reads shorter than 50 bp
   ```

2. Place input files in the appropriate subdirectories under `input/raw_reads`. The files should look like this:
   `{sample_name}_R1.fastq.gz`
   `{sample_name}_R2.fastq.gz`

3. Place already assembled genomes in the appropriate subdirectories under `input/genomes`. The files should look like this:
   `{genome_name}.fasta`

4. Set optional tool to true if you want to run it...

## Usage

### Running on the SIT Cluster using Slurm

```bash
module load anaconda3
module load mamba
conda activate snakemake
snakemake -s workflow/Snakefile.py --workflow-profile workflow/profiles/ga_pipeline/
```

### Running Locally

```bash
conda activate snakemake
snakemake -s workflow/Snakefile.py --workflow-profile workflow/profiles/default
```

Note: Remove the `-n` flag after verifying the dry run.

## Output

### Directory
The pipeline generates the following outputs:
```
.
├── 01_qc
│   ├── fastp
│   ├── fastqc
│   ├── multiqc
│   ├── quast
│   └── summary_qc.tsv
├── 02_kegganog
├── 03_pangenome
│   ├── pirate
│   └── genes_of_interest_presence_absence.tsv
├── 04_variant_calling
│   ├── snippy
│   ├── snippy-core
│   └── genes_of_interest_snp_matrix.tsv
├── 05_gwas
├── 06_typing
│   ├── chewbacca
│   ├── mlst
│   └── spaTyper
├── 07_amr
├── 08_temp
│   ├── temp_faa
│   ├── temp_fasta
│   ├── temp_gff
│   ├── temp_vcf
│   └── temp_vcf_gz
├── C22
│   ├── annotation
│   │   ├── bakta
│   │   ├── emapper
│   │   └── prokka
│   └── assembly
│       ├── contigs_fixed.fasta
│       ├── spades
│       └── unicycler
├── C24
├── C28
└── C32
```

## Some solutions

### Display anvi'o ringplot on your local computer
- Install anvi'o environment on your local computer.
- Download output/03_pangenome/anvio/anvio_Pangenome-PAN.db and output/03_pangenome/anvio/anvio_storage-GENOMES.db files.
- Display data:
```
conda activate anvio-8
anvi-display-pan -g anvio_storage-GENOMES.db -p anvio_Pangenome-PAN.db
```
-->make your ringplot pretty on the anvio server... (this could help: https://merenlab.org/2016/11/08/pangenomics-v2/)

### Download cgmlst alleles
- Download alleles as fasta from into `path/databases/{data_source_speciesname}/alleles/`, add in `config.yaml --> cgmlst: true` before running the pipeline.
	- cgmlst.org
	- chewbbaca.online
