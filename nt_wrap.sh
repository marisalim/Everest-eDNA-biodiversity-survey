#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N ntdb
#PBS -q long-24core

module load shared
module load blast+/2.7.1

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
gunzip -c nt.gz | makeblastdb -in - -dbtype nucl -out ntdb -title ntdbtitle -parse_seqids
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -xzvf taxdb.tar.gz
