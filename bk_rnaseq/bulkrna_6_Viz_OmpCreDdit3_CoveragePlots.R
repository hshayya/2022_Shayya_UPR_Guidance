#Create Normalized Bigwigs 

library(Rsubread)
library(rtracklayer)
library(tidyverse)
library(reshape2)
library(patchwork)
library(DESeq2)

###########################################################
###        Prep for Normalized Bigwig Generation        ###
###########################################################

#Will use --scaleFactor 1/size factor from DESeq2 for bigwig normalization in bamcoverage
#Ref: https://github.com/deeptools/deepTools/issues/401

#Compute size factors from DESeq2.

#Start w/ list of bam file locations. See STAR alignment code in "alig" folder. 
bam_files <- 
  c('/path/to/bams/batch/a',
    '/path/to/bams/batch/b',
    '/path/to/bams/batch/c') %>%
  map(list.files, pattern = '.*q30.bam$', recursive = T, full.names = T) %>%
  unlist() 

#Count Reads
cts_star_out <- map(bam_files, .f = function(x) {
  Rsubread::featureCounts(files = x,
                          annot.ext = 'IS_merged_CodingTranscriptsHS.gtf',
                          isGTFAnnotationFile = T,
                          GTF.featureType = 'exon',
                          allowMultiOverlap = T,
                          useMetaFeatures = T,
                          readShiftType = 'downstream',
                          readShiftSize = 0,
                          read2pos = NULL, #no read reduction, allow any part of read to overlap
                          strandSpecific = 2, 
                          isPairedEnd = T, 
                          nthreads = 15)$counts
})

#Create Matrix
cts_star_matrix <- map(cts_star_out, as.data.frame) %>% bind_cols()
colnames(cts_star_matrix) <- str_extract(bam_files, 'OmpCre_(WT|Ctrl|Exp)[0-9]+') #on my files I call WT=WT, Ctrl = Het, Exp = cKO. If aligned from GEO you'll already have the WT/Het/cKO names updated.
cts_star_matrix <- as.matrix(cts_star_matrix)

#Batch Data. 
#Again note WT/Ctrl/Exp (my internal terminology) = WT/Het/cKO on GEO
batch_mappings <- bind_rows(
  data.frame('lib_nm' = c('Exp2','Exp3','Exp4','Ctrl3'),
             'batch' = 'a'), #seq 10/22/21
  data.frame('lib_nm' = c('Exp1','Ctrl1','Ctrl2'),
             'batch' = 'b'), #seq 11/8/21 (submitted before a, but seq took forever- hence apparently backwards order!)
  data.frame('lib_nm' = c('Exp5','Ctrl4','Ctrl5','WT1','WT2','WT3','WT4'),
             'batch' = 'c')
)

#Coldata
star_coldata <- data.frame('lib_nm' = colnames(cts_star_matrix)) %>%
  left_join(mutate(batch_mappings, lib_nm = paste0('OmpCre_', lib_nm)), 
            by = 'lib_nm') %>%
  mutate('gt' = str_extract(lib_nm,'WT|Ctrl|Exp'),
         'batch' = fct_relevel(batch,'a','b'),
         'gt' = fct_relevel(gt, 'WT','Ctrl')) %>%
  column_to_rownames(var = 'lib_nm')
  
star_coldata <- coldata_ompcreonly[colnames(cts_star_matrix),]

dds_star <- DESeqDataSetFromMatrix(countData = cts_star_matrix,
                                   colData = star_coldata,
                                   design = ~batch + gt) %>%
  DESeq()

#Parameters for bamcoverage
#input bam file | output file path prefix | 1/sizefactor
bamcoverage_params <- data.frame(
  'input' = bam_files,
  'output' = map_chr(bam_files, .f = function(x) {
    base_ <- str_split(x, '/')[[1]]
    f_final <- paste0(str_extract(base_[length(base_)], '.*(?=.bam)'), '_norm_')
    #paste0(c(base_[-length(base_)], f_final), collapse = '/') #to put files in native folders
    paste0('/path/to/base/dir/for/norm/bw/files',
           f_final) #put in one designated folder!
  }),
  'norm_factor' = 1/sizeFactors(dds_star)
)

#write_tsv(bamcoverage_params, '/path/to/fout/for/parameter/file', col_names = F)

#At this point, use output tsv above to run bamCoverage (see bulkrna_6a_run_bamcoverage.sh code), read bw in below to continue with analysis

####################################################
###        Generate Ddit3 Coverage Figure        ###
####################################################

#Parse GTF of gene models and filter for Ddit3
gtf_ranges <- rtracklayer::import('IS_merged_CodingTranscriptsHS_sorted.gtf')

ddit3_granges <- gtf_ranges[gtf_ranges$transcript_id == 'ENSMUST00000026475']

#Found the loxP site locations by looking at Online Figure 1 from https://pubmed.ncbi.nlm.nih.gov/25872946/. 
#They have the sequence shown right by the loxP sites there, so I could just figure it out by looking at intron sequences in Ensembl
loxp_5 <- end(ddit3_granges[ddit3_granges$type == 'exon'])[2] + 238 #5'loxP site is ~238bp into intron2 (end of second exon)
loxp_3 <- end(ddit3_granges[ddit3_granges$type == 'exon'])[3] + 131 #3'loxP site is ~131bp into intron3 (end of third exon)

ddit3_genomic_range <- GRanges(seqnames = 'chr10',
                               ranges = IRanges(start = min(start(ddit3_granges))-100,
                                                end = max(end(ddit3_granges)) + 100))

#Parse Normalized Bigwigs
ompcre_ddit3_bwtraces <- 
  list.files('/path/to/base/dir/for/norm/bw/files',
             pattern = 'fwd.bw', full.names = T) %>%
  map_dfr(.f = function(x) {
    bw_dat <- rtracklayer::import(x,which = ddit3_genomic_range)
    lab_ <-  str_split(x,'/')[[1]] %>% .[length(.)] %>% str_extract('(WT|Ctrl|Exp)[0-9]+')
    
    as.data.frame(bw_dat@ranges) %>%
      mutate('y_min' = 0,
             'score' = bw_dat$score,
             'lab' = lab_,
             'gt' = str_extract(lab_, 'WT|Ctrl|Exp'))
  })

#Compute per-nt average across genotypes 
#Note this is fairly slow: bw has position encoded in blocks (ie. start, end , y val in that block). 
#Here we have to collapse to per-nt position, then average and plot 
avg_bwtraces <- 
  pmap_dfr(ompcre_ddit3_bwtraces, .f = function(start,end,width,y_min,score,lab,gt) {
    data.frame('x' = seq(start-1,end-1), #end-1 here b/c not using "end", this is bars
               'score' = score,
               'lab' = lab,
               'gt' = gt)
  }) %>%
  group_by(x, gt) %>%
  summarize(score = mean(score))

#Final Figure
{
  #Bigwig Traces
  ggplot(data = avg_bwtraces %>%
          mutate('gt' = fct_recode(gt, 'Het' = 'Ctrl','cKO' = 'Exp'),
                 'gt' = fct_relevel(gt, 'WT','Het')),
        aes(x = x, y = score,
            fill = gt)) +
    geom_area() +
    facet_wrap(facets =vars(gt), ncol = 1) +
    geom_vline(xintercept = loxp_5, linetype = 'dashed') +
    geom_vline(xintercept = loxp_3, linetype = 'dashed') +
    theme_bw() +
    scale_y_continuous(breaks = c(0,1500)) + 
    theme(axis.text.x = element_blank(), plot.background = element_blank(), 
          panel.grid = element_blank(), axis.title = element_blank(),
          axis.ticks.x = element_blank(), strip.background = element_blank(),
          strip.text= element_blank(), panel.border = element_blank()) +
    guides(fill = guide_legend(title = 'Geno'))} / 
    {
      #Annotation
      ggplot(data = as.data.frame(ddit3_genomic_range) %>%
               mutate('y_' = 1, 'y_end' = 1), 
             aes(x = start-1, xend = end, y = y_, yend = y_end)) +
        geom_blank() +
        geom_segment(data = as.data.frame(ddit3_granges[ddit3_granges$type == 'exon']@ranges),
                     y = 1, yend = 1, size = 2) +
        geom_segment(data = as.data.frame(ddit3_granges[ddit3_granges$type == 'CDS']@ranges),
                     y = 1, yend = 1, size = 5) +
        geom_segment(data = data.frame('start' = min(start(ddit3_granges[ddit3_granges$type == 'exon'])),
                                       'end' = max(end(ddit3_granges[ddit3_granges$type == 'exon']))),
                     y = 1, yend = 1, size = 1) +
        geom_point(data = data.frame(x = c(loxp_5, loxp_3)),
                   aes(x = x), color = 'red', y = 1, shape = 17, size = 3, inherit.aes = F) +
        theme_void()
    } + plot_layout(heights = c(0.8,0.2))

