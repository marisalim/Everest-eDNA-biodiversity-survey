#!/bin/bash
#SBATCH -p med2
#SBATCH -J 16s_r1bat1
#SBATCH --mem=10G
#SBATCH -t 30:00:00
#SBATCH --mail-type=ALL

source activate singlem

REF='/singlem_test/singlem_extra_packages/release1/4.40.2013_08_greengenes_97_otus.with_euks.spkg'
SEQ_PATH='./raw_data'
OUT_PATH='./results_16s_singlem'
READ='R1'

# R1 file, loop
for i in {13..15}; do echo ${i}; singlem pipe --singlem_packages ${REF} --forward ${SEQ_PATH}/${i}_${READ}_001.fastq.gz --output_extras --otu_table ${OUT_PATH}/ts${i}${READ}_16sotu_table.csv --threads 4; done
for i in {17..33}; do echo ${i}; singlem pipe --singlem_packages ${REF} --forward ${SEQ_PATH}/${i}_${READ}_001.fastq.gz --output_extras --otu_table ${OUT_PATH}/ts${i}${READ}_16sotu_table.csv --threads 4; done
