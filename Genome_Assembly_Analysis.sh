#!/bin/bash

# BINF 6110 Assignment 1 - Analysis Pipeline
# Analysis of Salmonella enterica genome assembly and comparison to reference genome

# Frist connected to compute canada 
# ssh mahnoorn@narval3.computecanada.ca
# Generate SSH key
# Add SSH key to GitHub
# Authenticate ssh -T git@github.com
# Clone repository git clone git@github.com:mahnoor-nizz/BINF6110-Assignment-1-.git


#Load Modules
module load StdEnv/2023
module load sra-toolkit/3.0.9
module load gcc/12.3
module load arrow/22.0.0
module load apptainer/1.4.5
module load minimap2/2.28
module load mummer/4.0.0
module load bcftools/1.22
module load samtools/1.22.1
module load python/3.11.5


# Set up directories
mkdir -p data qc assembly alignment variants visualization

# Download data
cd data
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565

wget "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000006945.2/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" -O ncbi_dataset.zip
unzip ncbi_dataset.zip
mv ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna reference.fasta


# Convert SRA to FASTQ
fasterq-dump SRR32410565
mv SRR32410565.fastq raw_reads.fastq
cd ..

# Python environment to use nanoplot and flye in compute canada
python -m venv bioinfo_env
source bioinfo_env/bin/activate

pip install --upgrade pip
pip install NanoPlot
pip install medaka
pip install biopython
pip install pandas
pip install matplotlib
pip install seaborn
pip install numpy
pip install quast
pip install flyey


# Quality Control
NanoPlot --fastq data/raw_reads.fastq --outdir qc --plots dot kde --threads 8


# Genome Assembly
flye --nano-hq data/raw_reads.fastq --out-dir assembly --threads 8 --genome-size 5m --iterations 2

# Assembly Quality Assessment
quast assembly/assembly.fasta -r data/reference.fasta -o qc/quast_results --threads 8

# Aligning Assembly to reference genome
# create index for reference genome
minimap2 -d data/reference.mmi data/reference.fasta
minimap2 -ax asm5 -t 8 --cs data/reference.mmi assembly/assembly.fasta > alignment/assembly_vs_ref.sam

# Align reads to reference genome for variant calling
minimap2 -ax map-ont -t 8 data/reference.fasta data/raw_reads.fastq > alignment/reads_vs_ref.sam
# Convert SAM to sorted BAM

# Convert Assembly alignment to BAM and sort
samtools view -bS alignment/assembly_vs_ref.sam > alignment/assembly_vs_ref.bam
samtools sort alignment/assembly_vs_ref.bam -o alignment/assembly_vs_ref.sorted.bam
samtools index alignment/assembly_vs_ref.sorted.bam
samtools flagstat alignment/assembly_vs_ref.sorted.bam > alignment/assembly_alignment_stats.txt
samtools coverage alignment/assembly_vs_ref.sorted.bam > alignment/assembly_coverage_summary.txt

# Convert read alignment to BAM and sort
samtools view -bS alignment/reads_vs_ref.sam > alignment/reads_vs_ref.bam
samtools sort alignment/reads_vs_ref.bam -o alignment/reads_vs_ref.sorted.bam
samtools index alignment/reads_vs_ref.sorted.bam
samtools flagstat alignment/reads_vs_ref.sorted.bam > alignment/reads_alignment_stats.txt
samtools depth alignment/reads_vs_ref.sorted.bam > alignment/reads_coverage.txt
samtools coverage alignment/reads_vs_ref.sorted.bam > alignment/reads_coverage_summary.txt

# index reference genome for variant calling
samtools faidx data/reference.fasta

# Variant Calling
bcftools mpileup -f data/reference.fasta alignment/reads_vs_ref.sorted.bam -Ou | \
bcftools call -mv --ploidy 1 -Ov -o variants/raw_variants.vcf

# Filter Variants
bcftools filter -i 'QUAL>=20 && DP>=10' variants/raw_variants.vcf -o variants/filtered_variants.vcf

# Stats
bcftools stats variants/filtered_variants.vcf > variants/variant_stats.txt

# structural Variant
dnadiff -p variants/comparison \
        data/reference.fasta \
        assembly/assembly.fasta
