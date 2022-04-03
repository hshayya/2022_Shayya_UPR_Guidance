library(tidyverse)
library(reshape2)
library(DESeq2)

#Load Files
setwd('/media/storageA/hani/Olfactory_Regulon')

#Commented-out code included for completeness only. Individual Quant.sf files not uploaded to Zenodo
#Just use the DEseq2 RDS object if running this.

#TxImport
#all_count_files <- list.files('/media/storageA/hani/Olfactory_Regulon', recursive = T, pattern = 'quant.sf', full.names = T)
#tx2gene <- read.table('/media/storageA/hani/Lhx2_UPR_Computation/tx_to_gene_name.tsv', sep='\t', 
#                      header=F, stringsAsFactors = F, col.names = c('tx','gene'))
#all_lineage_tximport <- tximport::tximport(files = all_count_files, type = 'salmon', tx2gene = tx2gene)

#Metadata
#colData_all_lineage <- tibble('file_names' = all_count_files) %>%
#  mutate('description' = map_chr(file_names, 
#                                 .f = function(x) {str_split(x, '/') %>% .[[1]] %>% .[length(.)-1] %>% str_extract('(?<=SalmonOut_).*')})) %>%
#  column_to_rownames(var = 'description')

#Create a DESeq2 Object From TxImport
#deseq_obj_all_lineage <- DESeq2::DESeqDataSetFromTximport(txi = all_lineage_tximport, 
#                                                          colData = colData_all_lineage, 
#                                                          design = ~1)
deseq_obj_all_lineage <- readRDS('272samples_Deseq2Object_RawCts.rds')


#Check for Libs w/ Alignment Problems
ggplot(data = colSums(counts(deseq_obj_all_lineage)) %>% log10() %>% {tibble('counts' = ., lib = names(.))}, 
       aes(x = counts)) +
  geom_line(stat = 'density') +
  theme_bw() +
  geom_point(aes(y = 1), alpha = 0.1)

#Remove libraries with <1000 counts and genes with >10 counts total across all libraries
deseq_obj_all_lineage <- deseq_obj_all_lineage[rowSums(counts(deseq_obj_all_lineage)) > 10,
                                               colSums(counts(deseq_obj_all_lineage))>10^3] 

#Run VST
vst_stabilized <- DESeq2::vst(deseq_obj_all_lineage)
#assays(vst_stabilized)[[1]] %>% as.data.frame() %>% rownames_to_column(var = 'gene') %>% 
#  write_tsv('./272samples_OE_expression_matrix_forArachne.tsv', col_names = T)


#Compute Tx Factor List for Aracne. Roughly inspired by https://www.ncbi.nlm.nih.gov/pubmed/27322546. 
ensembl <- biomaRt::useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', mirror = 'www')
ms_tx_factors <-biomaRt::getBM(attributes = c('external_gene_name', 'go_id','name_1006'),
                               values = c('GO:0003700','GO:0003677', 'GO:0006351', 'GO:0008134'), 
                               filters = 'go', mart = ensembl)

#Take either DNA-binding transcription factor activity or DNA bind and transcription, DNA-templated or transcription factor binding
ms_tx_list <- ms_tx_factors %>% filter(go_id %in% c('GO:0003700','GO:0003677', 'GO:0006351', 'GO:0008134')) %>%
  dcast(external_gene_name ~ go_id) %>%
  mutate_at(.vars = vars(-external_gene_name), .funs = function(x) {!is.na(x)}) %>%
  filter(`GO:0003700` | (`GO:0003677` & (`GO:0006351`|`GO:0008134`))) %>%
  pull(external_gene_name) %>% 
  .[!str_detect(.,'^(Hdac|Pol)')] #a couple of things I don't think are actually tfs
#ms_tx_list <- readLines('tf_list.txt')