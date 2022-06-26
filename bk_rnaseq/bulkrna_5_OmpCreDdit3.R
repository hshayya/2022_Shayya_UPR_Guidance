#Analyze the Omp-Cre Ddit3 Experiment

library(tidyverse)
library(reshape2)
library(DESeq2)
library(tximport)
library(patchwork)

#Imports
tx2gene <- read.table('tx_to_gene_name_gfp_tdtom_lacz.tsv', sep='\t', 
                      header=F, stringsAsFactors = F, col.names = c('tx','gene'))

guidance_network_toplot <- read_tsv('guidance_network.tsv')

iRFP_DE_High_vs_Low <- read_tsv('DE_iRFPBright_vs_Dim.tsv')

#Look at all the Omp-Cre Ddit3 Data
omp_cre_paths <- list.files(pattern = 'OmpCre.*sf') %>%
  set_names(str_extract(.,'.*(?=_quant.sf)')) 

imported_data_ompcreonly <- tximport::tximport(omp_cre_paths, type = 'salmon', tx2gene = tx2gene) 

batch_mappings <- bind_rows(
  data.frame('lib_nm' = c('cKO2','cKO3','cKO4','Het3'),
             'batch' = 'a'), #seq 10/22/21
  data.frame('lib_nm' = c('cKO1','Het1','Het2'),
             'batch' = 'b'), #seq 11/8/21 (submitted before a, but seq took forever- hence apparently backwards order!)
  data.frame('lib_nm' = c('cKO5','Het4','Het5','WT1','WT2','WT3','WT4'),
             'batch' = 'c')
)

coldata_ompcreonly <- data.frame('lib_id' = names(omp_cre_paths)) %>%
  mutate('gt' = str_extract(lib_id,'WT|Het|cKO'),
         'gt' = fct_relevel(gt, 'WT','Het'),
         'lib_nm' = str_extract(lib_id, '(WT|Het|cKO)[0-9]+')) %>%
  left_join(batch_mappings, by = 'lib_nm') %>%
  dplyr::select(lib_id, gt, batch) %>%
  column_to_rownames(var = 'lib_id')

dds_ompcreonly <- DESeq2::DESeqDataSetFromTximport(txi = imported_data_ompcreonly, 
                                                   colData = coldata_ompcreonly,
                                                   design = ~batch + gt) %>%
  DESeq2::DESeq(test = 'LRT', reduced = ~batch)

#PCA Plot
pca_ompcreonly_data <- as.data.frame(DESeq2::plotPCA(DESeq2::vst(dds_ompcreonly), 
                                                     intgroup = c('gt','batch'), returnData = T))

ggplot(data = pca_ompcreonly_data %>%
         mutate('gt' = fct_relevel(gt, 'WT','Het'), 
                'batch' = fct_relevel(batch,'a','b')),
       aes(x = PC1, y = PC2, color = batch, shape = gt)) +
  geom_point(size = 4) +
  xlab(paste0('PC1-', round(attr(pca_ompcreonly_data,'percentVar')[1]*100,0),'% Variance')) +
  ylab(paste0('PC2-', round(attr(pca_ompcreonly_data,'percentVar')[2]*100,0),'% Variance')) +     
  theme_bw() +
  labs(color = 'Batch', shape='Geno') +
  theme(legend.position = 'bottom') + 
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

#Extract results: lfcs cKO vs. Het, pvals by LRT.
res_ompcre_ddit3_lrt <- DESeq2::results(dds_ompcreonly, 
                                        contrast = list(c('gt_cKO_vs_WT'))) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  mutate('significance' = ifelse(padj < 0.05, 'Significant','Not Significant'))

#write_tsv(res_ompcre_ddit3_lrt, 'DE_OmpCreDdit3_cKOvsWT_LRT.tsv')

res_ompcre_ddit3_het_wt_lrt <- DESeq2::results(dds_ompcreonly, contrast = list(c('gt_Het_vs_WT'))) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  mutate('significance' = ifelse(padj < 0.05, 'Significant','Not Significant'))

#write_tsv(res_ompcre_ddit3_het_wt_lrt, 'DE_OmpCreDdit3_HetvsWT_LRT.tsv')

#MA Plot for Supplement, Omp-Cre Ddit3 cKO vs. WT
ggplot(data = res_ompcre_ddit3_lrt %>%
         mutate('identity' = ifelse(str_detect(gene_name, 'Olfr'), 
                                    'Olfactory Receptor','Other')) %>%
         filter(!is.na(padj)), #remove the independenly filtered genes (esp imp for limits calculation) 
       aes(x = log10(baseMean), y = log2FoldChange, color = significance)) +
  geom_point(aes(alpha = significance)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_wrap(facets = vars(identity), nrow = 1) +
  scale_color_manual(values = c('Significant' = 'red','Not Significant' = 'black')) +
  ggrepel::geom_label_repel(data = . %>% 
                              filter(gene_name %in% c('Kif5a','Cpn1') |
                                       (gene_name %in% guidance_network_toplot$gene_name & 
                                          significance == 'Significant')),
                            aes(label = gene_name), min.segment.length = 0, show.legend = F,
                            size = 3, box.padding = 0.5) +
  theme_bw() +
  scale_alpha_manual(values = c('Significant' = 0.6,'Not Significant' = 0.1)) +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(alpha = 1)), nrow = 1) +
  ylab('log2FC Omp-Cre\nDdit3 cKO vs. WT') + xlab('log10 Mean Normalized Counts')

#Volcano Plot, for Paper (using Ddit3 cKO vs. WT lfcs)
ggplot(data = res_ompcre_ddit3_lrt %>%
         left_join(dplyr::select(guidance_network_toplot, gene_name, stress_ident),
                   by = 'gene_name') %>%
         mutate('stress_ident' = replace_na(as.character(stress_ident),'Other')) %>%
         filter(!is.na(padj)) %>%
         mutate('stress_ident' = fct_relevel(stress_ident,'Other')) %>%
         arrange(stress_ident),
       aes(x = log2FoldChange, y = -log10(padj), color = stress_ident, alpha = stress_ident)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  scale_color_manual(values = c('High Stress' = 'red','Low Stress' = 'blue','Other' = 'black')) +
  scale_alpha_manual(values = c('Other' = 0.2,'High Stress' = 1, 'Low Stress' = 1), guide = 'none') +
  xlim(c(-1,1) * max(abs(res_ompcre_ddit3_lrt[!is.na(res_ompcre_ddit3_lrt$padj),]$log2FoldChange))) +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  guides('color' = guide_legend(nrow = 1, override.aes = list(alpha = 1))) + 
  xlab('log2FC Omp-Cre Ddit3 cKO vs. WT')

#Show log2FCs as a box plot
ompcreddit3_guidance_combined <- res_ompcre_ddit3_lrt %>%
  left_join(res_ompcre_ddit3_het_wt_lrt, by = c('gene_name','baseMean','stat','pvalue','padj','significance'),
            suffix = c('_cKOvsWT','_HetvsWT')) %>%
  right_join(dplyr::select(guidance_network_toplot, gene_name,
                           log2FoldChange_scRNA = log2FoldChange), by = 'gene_name') %>%
  dplyr::arrange(desc(log2FoldChange_scRNA)) %>%
  mutate('gene_name' = factor(gene_name, levels = gene_name),
         'significance' = replace_na(significance,'Not Significant')) %>%
  melt(id.vars = c('gene_name','significance'), 
       measure.vars = c('log2FoldChange_HetvsWT','log2FoldChange_cKOvsWT','log2FoldChange_scRNA')) %>%
  mutate('variable' = fct_recode(variable, 'OmpCre_HetvsWT' = 'log2FoldChange_HetvsWT',
                                 'OmpCre_cKOvsWT' = 'log2FoldChange_cKOvsWT',
                                 'scRNA_HighvsLow' = 'log2FoldChange_scRNA'),
         'variable' = fct_relevel(variable,'OmpCre_cKOvsWT','OmpCre_HetvsWT','scRNA_HighvsLow')) %>%
  dplyr::rename(log2FC = value) %>%
  dplyr::arrange(significance) #plot significant boxes last (so not overwritten by non-sig edges)

ggplot(data = ompcreddit3_guidance_combined %>%
         mutate_at(.vars = vars(gene_name,variable), .funs = fct_relevel, rev),
       aes(x = variable, y = gene_name, color = significance)) +
  geom_blank() + #force computation of scales w/o restricting to data subsets.
  geom_tile(data = . %>%
              filter(variable == 'scRNA_HighvsLow'),
            aes(fill = log2FC), size = 0.5, width = 1) +
  scale_fill_gradient2(low = scales::hue_pal()(2)[2], mid = 'lightgrey', high = scales::hue_pal()(2)[1],
                       midpoint = 0, breaks = c(-0.5,0,0.5)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = . %>%
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
  coord_equal()

#Point Plot Log2FC Omp-Cre Ddit3 vs. Log2FC scRNAseq
ggplot(data = res_ompcre_ddit3_lrt %>%
         right_join(guidance_network_toplot, by = 'gene_name', 
                    suffix = c('_ddit3','_scRNAseq')),
       aes(x = log2FoldChange_scRNAseq, y = log2FoldChange_ddit3)) +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_smooth(method = 'lm', se = T, color = 'black', alpha = 0.1) +
  geom_point(aes(color = stress_ident)) +
  ggrepel::geom_label_repel(data = . %>% 
                              mutate('lab' = ifelse(significance_ddit3 == 'Significant' &
                                                      log2FoldChange_scRNAseq < 0,
                                                    gene_name,'')),
                            aes(label = lab, color = stress_ident), 
                            min.segment.length = 0, show.legend = F, 
                            xlim = c(NA, 0), ylim = c(NA,0), size = 3, box.padding = 1) +
  ggrepel::geom_label_repel(data = . %>% 
                              mutate('lab' = ifelse(significance_ddit3 == 'Significant' &
                                                      log2FoldChange_scRNAseq > 0,
                                                    gene_name,'')),
                            aes(label = lab, color = stress_ident), 
                            min.segment.length = 0, show.legend = F,
                            xlim = c(0,NA), size = 3, box.padding = 1) +
  theme_bw() +
  scale_color_manual(values = c('High Stress' = 'red','Low Stress' = 'blue')) +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 1)) +
  xlab('log2FC High vs. Low Stress\nSingle Cells') + ylab('log2FC Omp-Cre Ddit3\ncKO vs. WT')

#Distribution of ORs in OmpCre-Ddit3 cKO vs. WT
ggplot(data = res_ompcre_ddit3_lrt %>%
         filter(str_detect(gene_name,'Olfr')) %>%
         left_join(dplyr::select(iRFP_DE_High_vs_Low,
                                 gene_name, 'log2FC_bright_dim' = 'log2FoldChange'),
                   by = 'gene_name') %>%
         mutate('iRFP_ID' = ifelse(log2FC_bright_dim > 0, 'High Stress', 'Low Stress')) %>%
         filter(!is.na(iRFP_ID)),
       aes(x = log2FoldChange)) +
  geom_line(aes(color = iRFP_ID, group = iRFP_ID), stat = 'density', alpha = 0.6) +
  geom_line(data = res_ompcre_ddit3_lrt %>% 
              filter(str_detect(gene_name, 'Olfr')),
            aes(color = 'All ORs'), stat  = 'density') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_bw() + theme(legend.title = element_blank(), legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  scale_color_manual(values = c('All ORs' = 'black',scales::hue_pal()(2) %>% set_names(c('High Stress','Low Stress')))) +
  xlab('log2FoldChange Omp-Cre Ddit3\ncKO vs. WT') + ylab('Density') +
  xlim(c(-4,4)) #cuts off one OR that has log2FC of 16.. 

#Show the Counts for individual replicates w/ the Axon Guidance Molecules
ggplot(data = res_ompcre_ddit3_lrt %>%
         filter(gene_name %in% guidance_network_toplot$gene_name & 
                  significance == 'Significant') %>%
         pull(gene_name) %>%
         {counts(dds_ompcreonly, normalized = T)[.,]} %>%
         as.data.frame() %>%
         rownames_to_column(var = 'gene_name') %>%
         melt(id.vars = c('gene_name'), variable.name = 'library') %>%
         left_join(rownames_to_column(as.data.frame(colData(dds_ompcreonly)), var = 'library'), 
                   by = 'library') %>%
         mutate('gene_name' = factor(gene_name, levels = dplyr::arrange(res_ompcre_ddit3_lrt, padj)$gene_name)),
       aes(x = gt, y = value, color = gt)) +
  geom_point(position = position_jitter(width = 0.2)) +
  facet_wrap(facets = vars(gene_name), nrow = 2, scales = 'free') +
  theme_bw() + ylab('Normalized Counts') + theme(legend.title = element_blank(),
                                                 axis.title.x = element_blank(), legend.position = 'none')
