#!/usr/bin/env python3
# ---------------------------------------------
# Marisa Lim
# Parse blast output for contigs and/or reads
# Extract query full sequence for top matches
# output table and plots
# Example bash loops:
# REF='chiro'
# ## for spades contigs
# for i in {13..33} 35; do python ./map_scripts/mapblast_parse.py --ref $REF --samp ts${i} --querytype contig --blastfile ./${REF}_map_results/ts${i}${REF}_blastout --queryseqfile ./${REF}_map_results/ts${i}_${REF}_spades/contigs.fasta; done
# ## for reads
# for i in {13..33} 35; do python ./map_scripts/mapblast_parse.py --ref $REF --samp ts${i} --querytype read --blastfile ./${REF}_map_results/ts${i}${REF}_reads_blastout --queryseqfile ./${REF}_map_results/ts${i}${REF}.fasta; done
#
# run:
# bash mapblast_batcher.sh
# March 2020
# ---------------------------------------------

import os, sys, argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description='''Parse blast output for contigs and/or reads
    ''',
    epilog='''Example: python mapblast_parse.py
    --ref chiro
    --samp ts13
    --querytype contig
    --blastfile ts13chiro_blastout
    # if spades contigs
    --queryseqfile contigs.fasta
    # if merged + unmerged reads
    --queryseqfile ts13chiro.fasta
    '''
)
parser.add_argument('--ref', help='ref species for BWA mapping', required=True)
parser.add_argument('--samp', help='sample ID', required=True)
parser.add_argument('--querytype', help='contig or read', required=True)
parser.add_argument('--blastfile', help='blast output file name', required=True)
parser.add_argument('--queryseqfile', help='contig or read input sequence file for blast search', required=True)

args=parser.parse_args()
arg_dict=vars(args)

## Filter blast output
if os.stat(arg_dict['blastfile']).st_size != 0:
    blastdat = pd.read_csv(arg_dict['blastfile'], sep='\t', header=None)
    blastdat.columns = ['queryID', 'genbankID', 'per_id', 'aln_len', 'num_mismatch', 'gap_open', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'qcovs', 'stitle']
    # print(blastdat.head())
    sampID = arg_dict['samp']
    # print(sampID)

    ## Subset table
    tophit_df = []
    uniq_blast_queries = blastdat.queryID.unique()
    # print(uniq_blast_queries)
    for myquery in uniq_blast_queries:
        # print(myquery)
        query_subset = blastdat.loc[blastdat['queryID'] == myquery].sort_values(by=['bitscore'], ascending=False)[0:1]
        stitle_len = len(query_subset.loc[:, 'stitle'].values[0].split(' '))
        print(stitle_len)
        if(stitle_len >= 3):
            query_subset['genus'] = query_subset.loc[:, 'stitle'].values[0].split(' ')[0]
            query_subset['species'] = query_subset.loc[:, 'stitle'].values[0].split(' ')[1]
            query_subset['subspecies'] = query_subset.loc[:, 'stitle'].values[0].split(' ')[2]
            print(query_subset)

            # Use query ID to grab sequence and add to table
            if os.stat(arg_dict['queryseqfile']).st_size != 0:
                seqdat = SeqIO.parse(open(str(arg_dict['queryseqfile'])), 'fasta')
                for fasta in seqdat:
                    name, sequence = fasta.id, str(fasta.seq)
                    # print(name) #only takes 1st element of header

                    if name == myquery:
                        query_subset['queryseq'] = sequence
            tophit_df.append(query_subset)
        elif(stitle_len <= 3):
            query_subset['genus'] = query_subset.loc[:, 'stitle'].values[0].split(' ')[0]
            query_subset['species'] = query_subset.loc[:, 'stitle'].values[0].split(' ')[0]
            query_subset['subspecies'] = query_subset.loc[:, 'stitle'].values[0].split(' ')[0]
            print(query_subset)

            # Use query ID to grab sequence and add to table
            if os.stat(arg_dict['queryseqfile']).st_size != 0:
                seqdat = SeqIO.parse(open(str(arg_dict['queryseqfile'])), 'fasta')
                for fasta in seqdat:
                    name, sequence = fasta.id, str(fasta.seq)
                    # print(name) #only takes 1st element of header

                    if name == myquery:
                        query_subset['queryseq'] = sequence
            tophit_df.append(query_subset)

    tophit_df2 = pd.concat(tophit_df)
    tophit_df2['ref'] = arg_dict['ref']
    tophit_df2['querytype'] = arg_dict['querytype']
    tophit_df2['sampleID'] = sampID

    ## output table
    tophit_df2.to_csv(str(arg_dict['blastfile'])+'_parsed.txt', sep='\t', index=False)

# if file is empty
elif os.stat(arg_dict['blastfile']).st_size == 0:
    print('Empty file. No blast results to parse.')
