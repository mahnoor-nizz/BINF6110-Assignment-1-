#!/bin/bash

# BINF 6110 Assignment 1 - Analysis Pipeline
# Analysis of Salmonella enterica genome assembly and comparison to reference genome

# Frist connected to compute canada 
# ssh mahnoorn@narval3.computecanada.ca
# Generate SSH key
# Add SSH key to GitHub
# Authenticate ssh -T git@github.com
# Clone repository git clone git@github.com:mahnoor-nizz/BINF6110-Assignment-1-.git


mkdir -p data qc assembly alignment variants visualization
cd data
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565

module load StdEnv/2023
module load sra-toolkit/3.0.9

fasterq-dump SRR32410565
mv SRR32410565.fastq raw_reads.fastq

wget "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000006945.2/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" -O ncbi_dataset.zip
unzip ncbi_dataset.zip
mv ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna reference.fasta
cd ..

# load python virtual environment for NanoPlot 
module load python/3.13
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


# Quality Control
NanoPlot --fastq data/raw_reads.fastq --outdir qc --plots dot kde --threads 4
