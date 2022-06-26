#Analyze the 4wk M28 Ddit3 Experiment

library(tidyverse)
library(reshape2)
library(DESeq2)
library(tximport)

#Imports
tx2gene <- read.table('tx_to_gene_name_gfp_tdtom_lacz.tsv', sep='\t', 
                      header=F, stringsAsFactors = F, col.names = c('tx','gene'))

guidance_network_toplot <- read_tsv('guidance_network.tsv')

###################################################################
################### 4wk M28 Ddit3 Full Analysis ###################
###################################################################

#Set up fpaths
m28_fpaths <- list.files(pattern = '4wkM28.*sf') %>%
  set_names(str_extract(.,'.*(?=_quant.sf$)'))

#Prep for Deseq2
m28_coldata <- data.frame('lib_id' = names(m28_fpaths)) %>%
  mutate('Genotype' = str_extract(lib_id, 'WT|cKO|Het'),
         'Genotype' = fct_relevel(Genotype, 'WT','Het'),
         'Population' = str_extract(lib_id, 'GFP|tdtom'),
         'Population' = fct_recode(Population,'tdtomato' = 'tdtom'),
         'Population' = fct_relevel(Population, 'GFP'),
         'ind' = str_extract(lib_id,'(?<=WT|Het|cKO)[0-9]+'),
         'ind' = fct_relevel(ind, function(x) {as.character(sort(as.numeric(x)))})) %>%
  column_to_rownames(var = 'lib_id')

all(names(m28_fpaths) == rownames(m28_coldata)) #looks good

m28_txi <- tximport::tximport(m28_fpaths, 
                              type = 'salmon', tx2gene = tx2gene) 

#Prep the Model Matrices for DESeq2
model_mat <- model.matrix(~ Genotype + Genotype:ind + Genotype:Population, m28_coldata)
model_mat <- model_mat[,apply(model_mat, 2, function(x) {!all(x==0)})]

model_mat_red <- model.matrix(~ Genotype  + Genotype:ind + Population, m28_coldata)
model_mat_red <- model_mat_red[,apply(model_mat_red, 2, function(x) {!all(x==0)})]

#Run DESeq2
dds_m28_lrt <- 
  DESeq2::DESeqDataSetFromTximport(txi = m28_txi, 
                                   colData = m28_coldata,
                                   design = ~1)

dds_m28_lrt@design <- model_mat

dds_m28_lrt <- dds_m28_lrt %>%
  DESeq2::DESeq(test = 'LRT', reduced = model_mat_red)

#PCA
pca_4wk_data <- DESeq2::plotPCA(DESeq2::vst(dds_m28_lrt), 
                                intgroup = c('Population','Genotype'), returnData = T)

ggplot(data = as.data.frame(pca_4wk_data) %>%
         mutate('Genotype' = fct_recode(Genotype, 'Het' = 'Ctrl','cKO' = 'Exp'),
                'Genotype' = fct_relevel(Genotype, 'WT'), 
                'Population' = fct_relevel(Population,'GFP')),
       aes(x = PC1, y = PC2, color = Population, shape = Genotype)) +
  geom_point(size = 3) +
  xlab(paste0('PC1-', round(attr(pca_4wk_data,'percentVar')[1]*100,0),'% Variance')) +
  ylab(paste0('PC2-', round(attr(pca_4wk_data,'percentVar')[2]*100,0),'% Variance')) +  
  scale_color_manual(values = c('GFP' = 'green','tdtomato' = 'red')) +
  guides(shape = guide_legend(title = 'Gt', nrow = 1),
         color = guide_legend(title = 'Population',nrow = 1)) +
  theme_bw() + theme(legend.position = 'bottom')

contrast_ = list(c('GenotypecKO.Populationtdtomato'), c('GenotypeWT.Populationtdtomato')) #KO vs. WT for model without interaction term
dds_m28_lrt_res <- DESeq2::results(dds_m28_lrt, contrast = contrast_) %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_name') %>%
  mutate('significance' = ifelse(padj < 0.05, 'Significant','Not Significant'))

#write_tsv(dds_m28_lrt_res, 'DE_M28_Ddit3_cKOvsWT_LRT.tsv')

#Volcano Plot KO vs WT using LRT p values!
ggplot(data = dds_m28_lrt_res %>% 
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
  xlab('Avg log2FoldChange tdt/GFP\nM28 Ddit3 cKO vs. WT') + ylab('-log10(padj)')

#MA Plot
ggplot(data = dds_m28_lrt_res %>% 
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
  ylab('Avg log2FoldChange tdt/GFP\nM28 Ddit3 cKO vs. WT') +
  xlab('log10 Mean Normalized Counts')

#Plot tdtom/gfp ratio for Nrp2 and Tenm2 Across Populations
ggplot(data = DESeq2::counts(dds_m28_lrt, normalized = T) %>%
         as.data.frame() %>%
         rownames_to_column(var = 'gene_name') %>%
         filter(gene_name %in% c('Nrp2','Tenm2')) %>%
         #filter(gene_name %in% c('Nrp2','Tenm4','Tenm2','Unc5b','Efna5','Robo2','Pcdh17','Epha7')) %>%
         melt(id.vars = c('gene_name'), variable.name = 'lib_id') %>%
         left_join(rownames_to_column(as.data.frame(colData(dds_m28_lrt)), var = 'lib_id'), 
                   by = 'lib_id') %>%
         group_by(gene_name, ind, Genotype) %>%
         summarize(value = log2(value[Population == 'tdtomato']/
                                  value[Population == 'GFP'])),
       aes(x = Genotype, y = value, color = Genotype)) +
  geom_point(position = position_jitter(width = 0.1)) +
  theme_bw() +
  facet_wrap(facets = vars(gene_name), scales = 'free') +
  theme(legend.position = 'none') +
  xlab('') + ylab('log2FoldChange Normalized Counts\ntdtomato/GFP')

#Compute the log2FCs for the het vs wt result and verify
res_het_vs_wt <-  DESeq2::results(dds_m28_lrt, contrast =  list(c('GenotypeHet.Populationtdtomato'), c('GenotypeWT.Populationtdtomato'))) %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_name') %>%
  mutate('significance' = ifelse(padj < 0.05, 'Significant','Not Significant'))

#write_tsv(res_het_vs_wt, 'DE_M28_Ddit3_HetvsWT_LRT.tsv')

#Show the Log2FCs as box plot
m28ddit3_guidance_combined <- dds_m28_lrt_res %>%
  left_join(res_het_vs_wt, by = c('gene_name','baseMean','stat','pvalue','padj','significance'),
            suffix = c('_cKOvsWT','_HetvsWT')) %>%
  right_join(dplyr::select(guidance_network_toplot, gene_name,
                           log2FoldChange_scRNA = log2FoldChange), by = 'gene_name') %>%
  dplyr::arrange(desc(log2FoldChange_scRNA)) %>%
  mutate('gene_name' = factor(gene_name, levels = gene_name)) %>%
  melt(id.vars = c('gene_name','significance'), 
       measure.vars = c('log2FoldChange_HetvsWT','log2FoldChange_cKOvsWT','log2FoldChange_scRNA')) %>%
  mutate('variable' = fct_recode(variable, 'M28Ddit3_HetvsWT' = 'log2FoldChange_HetvsWT',
                                 'M28Ddit3_cKOvsWT' = 'log2FoldChange_cKOvsWT',
                                 'scRNA_HighvsLow' = 'log2FoldChange_scRNA'),
         'variable' = fct_relevel(variable,'M28Ddit3_cKOvsWT','M28Ddit3_HetvsWT','scRNA_HighvsLow')) %>%
  dplyr::rename(log2FC = value) %>%
  dplyr::arrange(significance) #plot significant boxes last (so not overwritten by non-sig edges)

ggplot(data = m28ddit3_guidance_combined %>%
         filter(variable == 'scRNA_HighvsLow'),
       aes(x = gene_name, y = variable, color = significance)) +
  geom_tile(aes(fill = log2FC), size = 0.5, width = 1) +
  scale_fill_gradient2(low = scales::hue_pal()(2)[2], mid = 'lightgrey', high = scales::hue_pal()(2)[1],
                       midpoint = 0, breaks = c(-0.5,0,0.5)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = m28ddit3_guidance_combined %>% 
              filter(variable != 'scRNA_HighvsLow'),
            aes(fill = log2FC), size = 0.5, width = 1) +
  scale_fill_gradient2(low = 'purple',mid = 'lightgrey', high = 'green', midpoint = 0) +
  scale_color_manual(values = c('Significant' = 'black', 'Not Significant' = NA),
                     guide = 'none') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), plot.margin = margin(0,0,0,0,'pt')) +
  coord_equal() +
  theme(legend.position = 'bottom')
