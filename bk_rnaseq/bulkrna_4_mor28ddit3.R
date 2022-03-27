#Mor28-iCre/iGFP Ddit3 RNA-seq Analysis

library(tidyverse)
library(reshape2)
library(DESeq2)
library(tximport)

setwd('./')
tx2gene <- read.table('tx_to_gene_name_gfp_tdtom_lacz.tsv', sep='\t', 
                      header=F, stringsAsFactors = F, col.names = c('tx','gene'))
guidance_network_toplot <- read_tsv('guidance_network.tsv') #see scRNAseq analysis

#Load in files
m28_4wk_paths <- list.files('./',
                            pattern = 'quant.sf', recursive = T, full.names = T) %>%
  .[!str_detect(.,'nofluor')] %>%
  .[str_detect(.,'SalmonOut_4wkM28')] %>%
  set_names(map_chr(., .f = function(x) {x %>% str_split('/') %>% .[[1]] %>% .[length(.)-1] %>% 
      str_extract('(?<=SalmonOut_).*')}))

m28_txi <- tximport::tximport(m28_4wk_paths, type = 'salmon', tx2gene = tx2gene) 

m28_coldata <- data.frame('lib_id' = names(m28_4wk_paths)) %>%
  mutate('gt' = str_extract(lib_id,'WT|Het|cKO'),
         'population' = str_extract(lib_id,'tdtomato|tdtom|GFP'),
         'population' = fct_collapse(population, tdtomato = c('tdtomato','tdtom')),
         'population' = fct_relevel(population,'GFP','tdtomato'),
         'gt' = ifelse(population == 'GFP','WT',gt), #collapse so do NOT need interaction term. Use Ctrl if WT has been dropped obviously!
         'gt' = fct_relevel(gt, 'WT','Het','cKO')) %>%
  column_to_rownames(var = 'lib_id')

#LRT, M28-Ddit3 4wk Data
dds_m28_lrt <- 
  DESeq2::DESeqDataSetFromTximport(txi = m28_txi, 
                                   colData = m28_coldata,
                                   design = ~population + gt) %>%
  DESeq2::DESeq(test = 'LRT',reduced = ~population)

#PCA Plot
pca_4wk_data <- DESeq2::plotPCA(DESeq2::vst(dds_m28_lrt), intgroup = c('gt','population'), returnData = T)
ggplot(data = as.data.frame(pca_4wk_data) %>%
         mutate('gt' = fct_relevel(gt, 'WT'), 
                'population' = fct_relevel(population,'GFP')),
       aes(x = PC1, y = PC2, color = population, shape = gt)) +
  geom_point(size = 3) +
  xlab(paste0('PC1-', round(attr(pca_4wk_data,'percentVar')[1]*100,0),'% Variance')) +
  ylab(paste0('PC2-', round(attr(pca_4wk_data,'percentVar')[2]*100,0),'% Variance')) +  
  scale_color_manual(values = c('GFP' = 'green','tdtomato' = 'red')) +
  guides(shape = guide_legend(title = 'Effective Gt'),
         color = guide_legend(title = 'Population')) +
  theme_bw()

#DE Testing. Log2FC cKO vs. WT, p value from LRT.
contrast_ = list(c('gt_cKO_vs_WT')) 
res_m28_ko_vs_wt <- DESeq2::results(dds_m28_lrt, contrast = contrast_) %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_name') %>%
  mutate('significance' = ifelse(padj < 0.05, 'Significant','Not Significant'))
#write_tsv(dds_m28_lrt_res, 'DE_M28_Ddit3_LRT.tsv')

#Volcano Plot KO vs WT using LRT p values!
ggplot(data = res_m28_ko_vs_wt %>% 
         left_join(dplyr::select(guidance_network_toplot, gene_name, stress_ident), 
                   by = 'gene_name') %>%
         mutate('stress_ident' = replace_na(as.character(stress_ident), 'Other'),
                'stress_ident' = fct_relevel(stress_ident,'Other')) %>%
         arrange(stress_ident),
       aes(x = log2FoldChange, y = -log10(padj),
           color = stress_ident, alpha = stress_ident)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  theme_bw() +
  ggrepel::geom_label_repel(data = . %>% filter(significance == 'Significant' & stress_ident != 'Other'),
                            aes(label = gene_name), fill = 'white', show.legend = F, min.segment.length = 0,
                            alpha = 1) +
  coord_cartesian(xlim = c(-10,10)) +
  scale_alpha_manual(values = c(1,1,0.1) %>% set_names(c('High Stress','Low Stress','Other'))) +
  scale_color_manual(values = c('red','blue','black') %>% set_names(c('High Stress','Low Stress','Other'))) +
  theme(legend.position = 'bottom', legend.title = element_blank()) + 
  guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  xlab('log2FoldChange\nM28 Ddit3 cKO vs. WT') + ylab('-log10(padj)')

#MA Plot
ggplot(data = res_m28_ko_vs_wt %>% 
         mutate('or_identity' = ifelse(str_detect(gene_name,'Olfr'), 'Olfactory Receptor','Other')) %>%
         filter(!is.na(significance)) %>%
         mutate('significance' = fct_relevel(significance, 'Not Significant')) %>%
         dplyr::arrange(significance),
       aes(x = log10(baseMean), y = log2FoldChange,
           color = significance)) +
  geom_point(aes(alpha = significance)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_bw() + theme(legend.position = 'bottom', legend.title = element_blank()) +
  ggrepel::geom_label_repel(data = . %>% filter((gene_name %in% guidance_network_toplot$gene_name &
                                                   significance == 'Significant') |
                                                  gene_name %in% c('Nrp2','Tenm2','Olfr1507')),
                            aes(label = gene_name), fill = 'white', show.legend = F, min.segment.length = 0) +
  facet_wrap(facets = vars(or_identity), nrow = 1) +
  scale_color_manual(values = c('Significant' = 'red','Not Significant' = 'black')) +
  scale_alpha_manual(values = c('Significant' = 0.6,'Not Significant' = 0.1)) + 
  guides('color' = guide_legend(override.aes = list(alpha = 1), nrow = 1)) +
  ylab('log2FoldChange M28 Ddit3\ncKO vs. WT')

#Plot for Nrp2 and Tenm2 Across Populations
ggplot(data = DESeq2::counts(dds_m28_lrt, normalized = T) %>%
         as.data.frame() %>%
         rownames_to_column(var = 'gene_name') %>%
         filter(gene_name %in% c('Nrp2','Tenm2')) %>%
         #filter(gene_name %in% c('Nrp2','Tenm4','Tenm2','Unc5b','Efna5','Robo2','Pcdh17','Epha7')) %>%
         melt(id.vars = c('gene_name'), variable.name = 'lib_id') %>%
         left_join(rownames_to_column(as.data.frame(colData(dds_m28_lrt)), var = 'lib_id')) %>%
         mutate(gene_name = fct_relevel(gene_name,'Nrp2','Tenm4','Unc5b','Robo2','Tenm2'),
                'grp' = ifelse(population == 'GFP','WT.GFP',paste0(gt,'.tdt')),
                'grp' = fct_relevel(grp, 'WT.GFP','WT.tdt','Het.tdt')),
       aes(x = grp, y = value, color = grp)) +
  geom_point(position = position_jitter(width = 0.1)) +
  theme_bw() +
  facet_wrap(facets = vars(gene_name), scales = 'free', nrow = 2) + theme(legend.position = 'none') +
  xlab('') + ylab('Normalized Counts')
