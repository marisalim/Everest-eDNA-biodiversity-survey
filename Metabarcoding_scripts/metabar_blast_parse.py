#!/usr/bin/env python3
# ---------------------------------------------
# Marisa Lim
# Parse blast output for contigs and/or reads
# Extract query full sequence for top matches
# output table and plots
# Example for looping bash cmds: 
# for i in {13..33} 35; do python ../BWA_mitogenome_mapping/map_scripts/metabar_blast_parse.py --blastfile ./qcov90_perid70/ts${i}_barcontig_blastout --queryseqfile ts${i}_metabar_contigs.fasta; done
# for i in {13..33} 35; do python ../BWA_mitogenome_mapping/map_scripts/metabar_blast_parse.py --blastfile ./qcov80_perid70/ts${i}_barcontig_blastout --queryseqfile ts${i}_metabar_contigs.fasta; done
# for i in {13..33} 35; do python ../BWA_mitogenome_mapping/map_scripts/metabar_blast_parse.py --blastfile ts${i}_krakenreads_blastout --queryseqfile ts${i}kraken.fasta; done
# March 2020
# ---------------------------------------------

import os, sys, argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description='''Parse blast output for contigs and/or reads
    ''',
    epilog='''Example: python metabar_blast_parse.py
    --blastfile ts13chiro_blastout
    --queryseqfile contigs.fasta
    '''
)

parser.add_argument('--blastfile', help='blast output file name', required=True)
parser.add_argument('--queryseqfile', help='contig or read input sequence file for blast search', required=True)

args=parser.parse_args()
arg_dict=vars(args)

## Filter blast output
if os.stat(arg_dict['blastfile']).st_size != 0:
    blastdat = pd.read_csv(arg_dict['blastfile'], sep='\t', header=None)
    blastdat.columns = ['queryID', 'genbankID', 'per_id', 'aln_len', 'num_mismatch', 'gap_open', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'qcovs', 'stitle']
    # print(blastdat.head())
    # sampID = arg_dict['blastfile'].split('_')[0]
    sampID = arg_dict['blastfile'].split('_')[1].split('/')[1]
    # print(sampID)

    ## Subset table
    tophit_df = []
    uniq_blast_queries = blastdat.queryID.unique()
    # print(uniq_blast_queries)
    for myquery in uniq_blast_queries:
        # print(myquery)
        query_subset = blastdat.loc[blastdat['queryID'] == myquery].sort_values(by=['bitscore'], ascending=False)[0:1]
        query_subset['genus'] = query_subset.loc[:, 'stitle'].values[0].split(' ')[0]
        query_subset['species'] = query_subset.loc[:, 'stitle'].values[0].split(' ')[1]
        query_subset['subspecies'] = query_subset.loc[:, 'stitle'].values[0].split(' ')[2]
        # print(query_subset)

        ## Use query ID to grab sequence and add to table
        if os.stat(arg_dict['queryseqfile']).st_size != 0:
            seqdat = SeqIO.parse(open(str(arg_dict['queryseqfile'])), 'fasta')
            for fasta in seqdat:
                name, sequence = fasta.id, str(fasta.seq)
                # print(name) #only takes 1st element of header
                contigcov = fasta.description.split(' ')[1].replace(',','')

                if name == myquery:
                    query_subset['queryseq'] = sequence
                    query_subset['readspercontig'] = contigcov
        tophit_df.append(query_subset)

    tophit_df2 = pd.concat(tophit_df)
    tophit_df2['sampleID'] = sampID
    tophit_df3 = tophit_df2.loc[:, ['sampleID', 'queryID', 'readspercontig', 'genbankID', 'per_id', 'aln_len', 'num_mismatch', 'evalue', 'bitscore', 'qlen', 'qcovs', 'stitle', 'genus', 'species', 'subspecies', 'queryseq']] #don't need the start/end values
    # print(tophit_df3.head())

    ## output table
    tophit_df3.to_csv(str(arg_dict['blastfile'])+'_parsed.txt', sep='\t', index=False)

# if file is empty
elif os.stat(arg_dict['blastfile']).st_size == 0:
    print('Empty file. No blast results to parse.')
