#Align OR CDS's
from __future__ import division
import dill
from collections import defaultdict
from plastid import BAMGenomeArray, Transcript, VariableFivePrimeMapFactory, GenomeHash, GTF2_TranscriptAssembler
import numpy as np
import pandas as pd
import copy
import itertools
import re
from pathos.multiprocessing import ProcessingPool as Pool
from collections import Counter
import twobitreader as twobit

#Read in the data and the genome
transcripts=list(GTF2_TranscriptAssembler('IS_merged_CodingTranscriptsHS.gtf')) #see ../anot
genome = twobit.TwoBitFile('genome.2bit') #not included due to size. mm10, UCSC-style chrom numbers

#Transcript dictionary
tx_dictionary = {tx.attr['transcript_id']:tx for tx in transcripts}

#Build dictionary by transcript (gene_id: list(transcript objects))
gene_dict = defaultdict(list)
for transcript in transcripts:
    gene_dict[transcript.get_gene()].append(transcript)

#Reverse Gene Dict (transcript_id: gene_id)
tx_by_gene_dict = {tx.attr['transcript_id']: tx.attr['gene_id'] for tx in transcripts}

#Create Lookup Dictionary (common name: gene_id)
lookup_dict = {i.attr['gene_name']:i.attr['gene_id'] for i in transcripts}

#Parse Output from OMP-GFP
OMPGFP_SalmonOut = pd.read_csv('ompgfp_pairedend_salmonquant.sf', sep='\t', header=0, index_col = False) #paired end, deep-seq'd Omp-GFP sample
OMPGFP_SalmonOut['gene_id'] = [tx_by_gene_dict[tx] for tx in OMPGFP_SalmonOut.Name]
OMPGFP_groups = OMPGFP_SalmonOut.groupby(['gene_id'])

#Filter for OR genes
OR_gene_ids = [v for (k,v) in lookup_dict.iteritems() if re.match('Olfr',k) is not None]

#Iterate through these and return the transcript with the highest TPM
OR_gene_by_selected_tx = {}
for gene_id in OR_gene_ids:
    quant_tx = OMPGFP_groups.get_group(gene_id)
    highest_tpm = quant_tx.sort_values(by = ['TPM'], ascending=False)['Name'].iloc[0]
    OR_gene_by_selected_tx[gene_id] = tx_dictionary[highest_tpm]

#Return the CDS for each transcript as a FASTA file
with open('OR_CDS.fa','w') as fout:
    for tx in OR_gene_by_selected_tx.itervalues():
        fout.write('>' + tx.attr['gene_name'] + '_' + tx.attr['gene_id'] + '_' + tx.attr['transcript_id'] + '\n')
        fout.write(tx.get_cds().get_sequence(genome).upper() + '\n')
    #
