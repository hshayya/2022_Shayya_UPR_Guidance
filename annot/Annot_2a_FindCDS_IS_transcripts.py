from plastid import GTF2_TranscriptAssembler
import numpy as np

#Read in the GTF Files
ensembl_gtf = list(GTF2_TranscriptAssembler('olfrs_ensembl.gtf'))
IS_gtf = list(GTF2_TranscriptAssembler('olfrs_IS.gtf'))

def annotate_CDS(unannotated_segment_chain, list_annotated_segmentchains):
    #get genename for our unannotated segment chain
    genename = unannotated_segment_chain.get_gene()
    
    #Generate a list of Gene Names for annotated segmentchains
    ensembl_names = []
    for x in list_annotated_segmentchains:
        ensembl_names.append(x.get_gene())
    
    #find indices of annotated segment chains corresponding to our unnanotated gene
    indices = []
    for i, elem in enumerate(ensembl_names):
        if (genename == elem):
            indices.append(i)
        
    
    #give up if we don't have any annotated segment chains for our unnanotated gene
    if (len(indices) == 0):
        return None
    
    #subset our list of annotated segment chains to just give those for the gene of interest 
    ensembl_forgivengene = [list_annotated_segmentchains[i] for i in indices]
    
    #Get CDS annotations for Ensembl Transcripts corresponding to gene of interest
    ensembl_cds_annotation = []
    for x in ensembl_forgivengene:
        ensembl_cds_annotation.append(x.get_cds())
    
    #give up if we don't have any annotated CDS's for our unannotated gene
    if (len(ensembl_cds_annotation) == 0):
        return None
    
    #Test if any of the identified CDS for gene of interest are fully contained within the CUFFORT transcript of interest. In theory this adjusts for splicing because it looks at all genomic positions, not just the end boundaries.
    covered_cds = [i for i in ensembl_cds_annotation if unannotated_segment_chain.covers(i)]
    
    #give up if none of the identified CDS are found in our unnanotated transcript
    if(len(covered_cds) == 0):
        return None
    
    lengths_cds = []
    for x in covered_cds:
        lengths_cds.append(x.get_length())
    
    largest_cds = covered_cds[np.argmax(lengths_cds)]
    
    #Otherwise, define maximum CDS bounds on our unannotated (now annotated!) transcript in GENOMIC coordinates. 
    starts = largest_cds.cds_genome_start
    stops = largest_cds.cds_genome_end
    
    #convert to segmentchain coordinates. Note that strandedness doesn't matter for minimum start and stop (above) because those are defined as genomic coordinates. However, we need to be careful about converting start/stop genomic coordinates to segmentchain coordinates. This depends on strand, as below
    
    chrom = unannotated_segment_chain.chrom
    strand = unannotated_segment_chain.strand #note that covered_cds already has ensured CDS is on same chr and strand as unnanotated transcript
    
    if (strand == '+'):
        start_segmentchain = unannotated_segment_chain.get_segmentchain_coordinate(chrom, starts, strand)
        stop_segmentchain = unannotated_segment_chain.get_segmentchain_coordinate(chrom, stops-1,strand)
        return unannotated_segment_chain.get_subchain(start_segmentchain, stop_segmentchain+1)
    
    if (strand == '-'): #must add 1 here becuase of how indexing works when comign from other side. We are 0-based, half open and I think the problem is the half open part, which causes the left-most genomic stop point to be 1 higher than it should be- screwing up the start for the transcript-based indexing.
        start_segmentchain = unannotated_segment_chain.get_segmentchain_coordinate(chrom, stops-1, strand)
        stop_segmentchain = unannotated_segment_chain.get_segmentchain_coordinate(chrom, starts,strand)
        return unannotated_segment_chain.get_subchain(start_segmentchain, stop_segmentchain+1)
    

segments = []
for i in IS_gtf:
    segments.append(annotate_CDS(i,ensembl_gtf))

#Transcripts that could not be annotated 
segments_problems = []
for i,text in enumerate(segments):
    try:
        if text == 'Error_CheckBounds' or text == None: 
            segments_problems.append(i)
        
    
    except TypeError:
        pass
    

annotated_segments = [i for i in segments if i is not None]

outputfile = open('IS_ORseqs_CDSannotated_new.tsv', 'w')
for i in annotated_segments:
    line = i.as_gtf(feature_type = 'CDS')
    outputfile.write(line)

outputfile.close()
