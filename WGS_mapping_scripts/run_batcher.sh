#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N bwa
#PBS -q long-24core

# need
module load shared
module load anaconda/2
module load blast+/2.7.1

REF='/gpfs/scratch/mclim/EverestMetaGenomics/BWA_MAPPING/mito_refs/Chironomus_tepperi_mitogenome.fasta'
REFNAME='chiro'

FQPATH='/gpfs/scratch/mclim/EverestMetaGenomics/METAGEN'
NCBI='/gpfs/scratch/mclim/EverestMetaGenomics/NCBI_blast_nt'

for i in {13..33} 35
do
  mkdir ts${i}${REFNAME}_files
  cd ts${i}${REFNAME}_files
  
  echo 'map reads...'
  echo '______________________'
  bash ../bwa_batcher.sh $REF $FQPATH/${i}_R1_001.fastq.gz $FQPATH/${i}_R2_001.fastq.gz TS${i}${REFNAME}
  
  echo 'filter sam file...'
  echo '______________________'
  awk '$4 != "0"' bwa_TS${i}${REFNAME}.sam > ts${i}_${REFNAME}_hits.sam
  
  echo 'convert sam to fastq...'
  echo '______________________'
  samtools fastq ts${i}_${REFNAME}_hits.sam > ts${i}${REFNAME}hits.fastq
  
  echo 'make spades assembly...'
  echo '______________________'
  spades.py -o ts${i}_${REFNAME}_spades --12 ts${i}${REFNAME}hits.fastq -t 1
  
  echo 'blastn nt search on spades contigs...'
  echo '______________________'
  blastn -db $NCBI/ntdb -max_target_seqs 10 -outfmt "6 std qlen qcovs stitle" -out ts${i}${REFNAME}_blastout -qcov_hsp_perc 90 -perc_identity 80 -query ./ts${i}_${REFNAME}_spades/contigs.fasta
  
  echo 'blastn nt search on filtered reads -- run when there are no Spades contigs...'
  echo '______________________'
  
  echo 'first, we will merge any paired reads to blastn them as one query...'
  bbmerge.sh in=ts${i}${REFNAME}hits.fastq out=ts${i}${REFNAME}_merged.fq outu=ts${i}${REFNAME}_unmerged.fq
  cat ts${i}${REFNAME}_merged.fq ts${i}${REFNAME}_unmerged.fq > ts${i}${REFNAME}_bbreads.fq

  cat ts${i}${REFNAME}_bbreads.fq | paste - - - - | sed 's/^@/>/g' | cut -f1-2 | tr '\t' '\n' | cut -d ' ' -f1 > ts${i}${REFNAME}.fasta
  blastn -db $NCBI/ntdb -max_target_seqs 10 -outfmt "6 std qlen qcovs stitle" -out ts${i}${REFNAME}_reads_blastout -qcov_hsp_perc 90 -perc_identity 80 -query ./ts${i}${REFNAME}.fasta

  echo 'parse blast results...'
  python ../mapblast_parse.py --ref ${REFNAME} --samp ts${i} --blastfile ts${i}${REFNAME}_blastout --queryseqfile ./ts${i}_${REFNAME}_spades/contigs.fasta --querytype contig
  python ../mapblast_parse.py --ref ${REFNAME} --samp ts${i} --blastfile ts${i}${REFNAME}_reads_blastout --queryseqfile ./ts${i}${REFNAME}.fasta --querytype read

  echo 'clean up and compress large files...'
  echo '______________________'
  gzip bwa_TS${i}${REFNAME}.sam
  
## only if bam or sorted files exist from bwa_batcher.sh  
  gzip bwa_TS${i}${REFNAME}.bam
  gzip bwa_TS${i}${REFNAME}.sorted
  cd ../
done






