#!/bin/bash 
#SBATCH -p med2
#SBATCH -J singlem17_21
#SBATCH --mem=5G
#SBATCH -t 30:00:00
#SBATCH --mail-type=ALL

source activate singlem

# change i values, do 3-5 samples at a time spread over a few batch scripts
for i in {17..21}; do singlem pipe --forward ./raw_data/${i}_R1_001.fastq.gz --reverse ./raw_data/${i}_R2_001.fastq.gz --otu_table ./results_singlem/ts${i}_otu_table.csv --output_extras --threads 4; done
