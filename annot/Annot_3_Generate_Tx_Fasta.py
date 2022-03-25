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

#Script
#Read in the data and the genome
transcripts=list(GTF2_TranscriptAssembler('IS_merged_CodingTranscriptsHS.gtf'))
genome = twobit.TwoBitFile('genome.2bit') #not included due to size. mm10 genome, UCSC-style chrom numbers

#Build dictionary by transcript
gene_dict = defaultdict(list)
for transcript in transcripts:
    gene_dict[transcript.get_gene()].append(transcript)

#Create Lookup Dictionary
lookup_dict = {i.attr['gene_name']:i.attr['gene_id'] for i in transcripts}

#Generate the Fasta index
fout = open('IS_merged_CodingTranscriptsHS_transcripts.fasta','w')
for tx in transcripts:
    fout.write('>' + tx.attr['transcript_id'] + '\n')
    fout.write(tx.get_sequence(genome) + '\n')

fout.close()
