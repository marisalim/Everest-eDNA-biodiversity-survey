# Everest eDNA Biodiversity Survey Analysis

This repository holds the scripts we used to analyze Everest WGS (for bacterial and eukaryote hits) and metabarcoding data. 

<!---These scripts were developed for and accompany the *in prep* manuscript:
*Estimating Biodiversity Using Environmental DNA Analysis Across the Tree of Life from Mount Everest’s Southern Flank*
Marisa C.W. Lim, Batya Nightingale, Charles Xu, Adam Solon, Nick Dragone, Steve Schmidt, Stephan Halloy, Alex Tait, Sandra Elvin, Aurora Elmore, Anton Seimon, and Tracie A. Seimon--->

## Contents
1. [Required software](#software)
1. [Scripts](#scripts)
1. [WGS eukaryote pipeline details](#wgspipeline)
1. [Metabarcoding blast script details](#metabarpipeline)
1. [Run WGS pipeline](#runwgs)
1. [Run metabarcoding blast script](#runmetabar)
1. [Compare datasets](#comparedats)

## Required software <a name="software"></a>

    bwa [Version: 0.7.15-r1140]
    samtools [Version: 1.9 (using htslib 1.9)]
    spades [SPAdes genome assembler v3.10.1]
    BBMerge [v38.22] (part of BBtools suite)
    blastn [Nucleotide-Nucleotide BLAST 2.7.1+]
    singleM [v0.13.2]

## Scripts <a name="scripts"></a>

These files were run on an HPC to map WGS reads to reference mitochondrial and chloroplast genomes. All python scripts written for Python3. 

**Disclaimer**: these are the scripts we used; they were developed on a MacOS and are not meant to be an out-of-the-box pipeline. If you would like to use them and have any questions, please post an [issue](https://github.com/marisalim/Everest-eDNA-biodiversity-survey/issues/new).

[WGS-bacteria](./SingleM_scripts):
- batch scripts for running default singleM ribosomal panel: `singlem_batch.sh`
- batch scripts for running singleM with Greengenes database: `singlem_16s_R1s_batch.sh` and `singlem_16s_R2s_batch.sh`

[WGS-eukaryotes](./WGS_mapping_scripts):
- BWA mapping commands are in `bwa_batcher.sh`
- Blast parsing commands are in `mapblast_parse.py`
- primary batch script is `run_batcher.sh`
- add taxonomic info, plot results: `wgs_blast_dataviz.r`
- an example [snakemake workflow for this pipeline](./WGS_mapping_scripts/Snakemake_ver)

[metabarcoding](./Metabarcoding_scripts/):
- Blast parsing commands are in `metabar_blast_parse.py`
- add taxonomic info, plot results: `metabar_blast_dataviz.r`

[Comparing results from metabarcoding vs. WGS](./Scripts_to_compare_dats):
- compare datasets: `compare_eDNAdat_splists.r`
- barplots and heatmaps: `heatmap_plts.r`

[Asymptotic regression model script](./everest.R)

## WGS eukaryote pipeline details <a name="wgspipeline"></a>

This pipeline was designed to specifically pull out Eukaryotic taxa from the WGS sequencing data. Due to the incomplete status of Eukaryotic taxonomic representation in reference sequence databases and from the Everest region, this is an imperfect approach. However, it works reasonably well for taxon discovery and exploration. It would require further tuning to improve completeness in taxon identification. Here are some details about what the various pipeline scripts are doing:

*1. Map to reference*
- map to mitochondrial reference or whatever is available (download from NCBI). The goal of this approach is to find Eukaryotic taxa. Since they represent a minority of the WGS data, the goal of this mapping step is to narrow the data down to reads more likely to match Eukaryote references of choice and filter out non-target reads (e.g., bacteria and other microorganisms).
  - note: this entire approach works better with small genomes like mitochondrial or chloroplast references. Quite slow on whole genomes.
  - e.g., mitochondrial refs take about 2 days or a little more to run all 22 samples
  - e.g., *Juniperus* chloroplast genome was better for pulling out plant reads, also had bacteria, so filtered blast results to only keep chloroplast gene hits

- map paired-end reads with `bwa mem`
- [`bwa_batcher.sh`](./WGS_mapping_scripts/bwa_batcher.sh) has the bwa commands
- optional step: (not run in pipeline code, but you can run this separately to view mappings in terminal)
    > convert .sam to .bam with `samtools faidx` and `samtools import`
    > sort BAM file with `samtools sort` and `samtools index`
    > view bam file in tview with `samtools tview`

*2. Subset SAM files for 'good' hits*

- `awk '$4 != "0"' [sample].sam > [sample]_hits.txt` to check if there are any reads that mapped, plus associated mapping score info. if the file is blank, then no mapping. Thus, I'm filtering on this field, and removing the rows that have 0's.

The $4 field is the 1-based leftmost mapping position of SAM files. Note the following:
>> "The first base in a reference sequence has coordinate 1. POS is set as 0 for an unmapped read without coordinate.  If POS is 0, no assumptions can be made about RNAME and CIGAR."

More about SAM files: 
- https://samtools.github.io/hts-specs/SAMv1.pdf
- https://davetang.org/wiki/tiki-index.php?page=SAM

*3. Convert filtered SAM file to interleaved fasta file*

use: `samtools fastq`

*4a. Assemble interleaved reads with spades*
- of the reads that mapped to reference, try to generate contigs for blast search
- this is running pretty quickly because there are usually not very many reads to assemble after the mapping filter
- doesn't work if there are too few reads to assemble, e.g., 2 reads
- see explanation of Spades and k-mer assembly: 
  - http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2019/01/genomics_tutorial_2019.pdf (e.g., pg 70-)
  - https://groups.google.com/forum/#!topic/abyss-users/47xOSCy9uQc

*4b. Option: merge reads if no spades contigs*
- use BBmerge to produce merged and unmerged reads
- for our 150bp paired-end WGS dataset, means you can get a max. 300bp merged read, but likely shorter, but better than the 150bp unpaired read

*5. Nucleotide blast search for spades contigs and/or merged & unmerged reads*
- **Edit per your study goals:** note different blast thresholds for percent sequence similarity (`-qcov_hsp_perc`) and percent query coverage (`-perc_identity`), output has top 10 hits per query (`-max_target_seqs`) on lines 40 and 50.

*6. Parse blast output per sample and reference, summarize*
- parse to keep best hit per query (sorted by bit score - mostly works, but in some cases, e.g., a tie, you need to manually check)

## Metabarcoding blast script details <a name="metabarpipeline"></a>
The metabarcoding data was primarily analyzed in Geneious to generate contigs. After blast search to NCBI nt database, the blast parsing step was conducted with [`metabar_blast_parse.py`](./Metabarcoding_scripts/metabar_blast_parse.py). Same parsing logic as for WGS data. 

## Running the WGS pipeline <a name="runwgs"></a>
1. Set up blast reference database. For our Everest paper, I set up a local installation of the nt database (unstable connection errors on HPC made remote run option erratic) with [this script](./nt_wrap.sh). Blast db set up takes a long time (many hours) because nt db is very large; make sure you have enough space to save the database files. Do the following:
```
mkdir NCBI_blast_nt
cd NCBI_blast_nt
# download of nt.gz takes about 40 minutes (compressed file is 65Gb)
# makeblastdb part took about an hr (db files (nhr, nin, nog, nsd, nsi, nsq; index nal) total = 73 Gb)
  
sbatch nt_wrap.sh
```

2. Edit input paths in the `run_batcher.sh` script. These 4 variables must be edited with correct path/name:
```
# path to reference genome fasta file
REF='EverestMetaGenomics/BWA_MAPPING/mito_refs/Chironomus_tepperi_mitogenome.fasta'
# abbreviated name for reference (used for output file naming)
REFNAME='chiro'
# path to sequence read fastq files
FQPATH='EverestMetaGenomics/METAGEN'
# path to NCBI blast database
NCBI='EverestMetaGenomics/NCBI_blast_nt'
```
3. Run script. The blast parsing step, in particular, takes a while, so I recommend running this script on an HPC.
- **Note on timing:** the entire pipeline can take 20-48+ hours to run with mitochondrial or chloroplast genome refs; takes upwards of a week or more for full genomes - the rate limiting step is my blast parser. This is not particularly efficient for whole genome refs (too many hits, so parsing is slow) but works quickly for small genomes like mito or chloroplast!

`sbatch run_batcher.sh` 

## Run metabarcoding blast script <a name="runmetabar"></a>

Ran [`metabar_blast_parse.py`](./Metabarcoding_scripts/metabar_blast_parse.py) as a bash loop, where `--blastfile` is the blast output (format 6) and the `--queryseqfile` is the contig or read input sequence file for blast search (used to add the sequence to parsed output file), e.g.,:
```
for i in {13..33} 35; do python metabar_blast_parse.py --blastfile ts${i}_barcontig_blastout --queryseqfile ts${i}_metabar_contigs.fasta; done
```

## Compare datasets <a name="comparedats"><a/>
At this point in the analysis, I did post-processing of the blast results in R to:
  - add additional taxonomy information using R package `taxize` to [WGS data](./WGS_mapping_scripts/wgs_blast_dataviz.r) and to [metabarcoding data](./Metabarcoding_scripts/metabar_blast_dataviz.r)
  - [compare Everest dataset results](./Scripts_to_compare_dats/compare_eDNAdat_splists.r)
  - [make heat maps to visualize OTU matrices](./Scripts_to_compare_dats/heatmap_plts.r)
