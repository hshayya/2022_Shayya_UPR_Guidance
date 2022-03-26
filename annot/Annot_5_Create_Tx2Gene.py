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

transcripts=list(GTF2_TranscriptAssembler('IS_merged_CodingTranscriptsHS.gtf'))
lookup_dict = {i.attr['transcript_id']:i.attr['gene_name'] for i in transcripts}

#No fluors
with open('tx_to_gene_name.tsv','w') as fout:
    for k,v in lookup_dict.iteritems():
        fout.write(k+'\t'+v+'\n')

#With Fluors
with open('tx_to_gene_name_gfp_tdtom_lacz.tsv','w') as fout:
    for k,v in lookup_dict.iteritems():
        fout.write(k+'\t'+v+'\n')
#add the manual mappings for fluorophores
with open('tx_to_gene_name_gfp_tdtom_lacz.tsv','a') as fout:
    fout.write('egfp\tegfp\ntdtomato\ttdtomato\nlacZ\tlacZ')
