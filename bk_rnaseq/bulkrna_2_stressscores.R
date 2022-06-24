#Compute Stress Scores from Bulk iRFP Bright/Dim Sorts
library(tidyverse)
library(reshape2)
library(DESeq2)
library(tximport)

#Necessary Files
zonal_annotation <- read_tsv('ORs-by-zone.txt', col_names = c('gene_name','zone'))
atf5ko_vs_wt <- read_csv('Atf5KO_CebpgKO_OmpGFP.csv')
matsunami_ORs <- c('oOR_Matsunami.csv','uOR_Matsunami.csv') %>%
  map(.f = function(path_) {
    name_ <- str_extract(path_, '(u|o)OR')
    read_csv(path_) %>% 
      pull(id) %>% list() %>% set_names(name_)
  }) %>% unlist(recursive = F)
matsunami_de <- read_tsv('matsunami_rtp_dko_vs_WT.tsv')
tx2gene <- read.table('tx_to_gene_name_gfp_tdtom_lacz.tsv', sep='\t', 
                      header=F, stringsAsFactors = F, col.names = c('tx','gene'))

#Prep for Deseq
file_paths <- list.files(path = './', 
                         pattern = '.*quant.sf', recursive = T) %>% 
  .[str_detect(.,'Ctrl[1-4]')] %>%
  set_names(str_extract(.,'.*(?=_quant.sf$)'))

imported_data <- tximport(file_paths, type = 'salmon', tx2gene = tx2gene)

coldata <- data.frame('lib_id' = names(file_paths)) %>%
  mutate('Population' = str_extract(lib_id,'iRFP(GFP)?'),
         'Intensity' = str_extract(lib_id,'Low|High')) %>%
  column_to_rownames(var = 'lib_id')

#Stress score = bright dim term from lm ~ population + intensity in CONTROL libs only (NOT swap)
dds <- DESeqDataSetFromTximport(txi = imported_data, 
                                colData = coldata,
                                design = ~Population + Intensity) %>%
  DESeq()

res <- results(dds, contrast = c('Intensity','High','Low')) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  mutate('significance' = ifelse(padj < 0.05, 'Significant','Not Significant'))

#Write data frames for visualization 
#write_tsv(res,'DE_iRFP_Bright_vs_Dim.tsv')
#write_tsv(res %>% filter(str_detect(gene_name,'Olfr')) %>% 
#            dplyr::select(gene_name, 'stress_score' = log2FoldChange),
#          'stress_scores.tsv')

###################################################
#### Generate the Stress Score Dot Plot Figure ####
###################################################

#Compute the KDE
kde <- res %>%
  filter(str_detect(gene_name, 'Olfr') & !is.na(log2FoldChange))  %>%
  {density(.$log2FoldChange)}

#Generate the Plot
new_color_scale_2_17 <- scales::hue_pal()(7)[c(1,3,5,6,7)] %>% set_names(c('M71','Class I','Mor23','P2','Mor28')) %>%
  c('other' = 'black')
new_alpha_scale_2_17 <- ifelse(names(new_color_scale_2_17) == 'other', 0.1, 1) %>% set_names(names(new_color_scale_2_17))

set.seed(1)
ggplot(data = res %>%
         filter(str_detect(gene_name, 'Olfr') & !is.na(log2FoldChange)) %>%
         mutate('kde_val' = map_dbl(log2FoldChange, .f = function(x) {
           kde_est <- kde$y[which.min(abs(kde$x - x))]
           0 + runif(1,-kde_est, kde_est)
         })) %>%
         left_join(zonal_annotation, by = 'gene_name') %>%
         left_join(data.frame('gene_name' = c('Olfr151','Olfr16','Olfr17', 'Olfr1507','Olfr545'),
                              'OR_name' = c('M71','Mor23','P2','Mor28','Class I'), stringsAsFactors = F)) %>%
         mutate('OR_name' = replace_na(OR_name, 'other'),
                'zone' = replace_na(as.character(zone), 'other'),
                'OR_name' = ifelse(OR_name == 'other', ifelse(zone == 'fishOR','Class I','other'), OR_name),
                'OR_name' = fct_relevel(OR_name,'Class I','M71','Mor23','P2','Mor28','other'),
                'OR_name' = fct_relevel(OR_name, rev)), 
       aes(x = kde_val, y = log2FoldChange, color = OR_name)) +
  geom_polygon(data = data.frame('log2FoldChange' = kde$x, 'kde_val' = kde$y) %>%
                 bind_rows(mutate(., 'kde_val' = -kde_val) %>% arrange(-log2FoldChange)), 
               fill = NA, color = 'black', alpha = 1) +
  geom_point(data = . %>% filter(OR_name == 'other'), aes(alpha = OR_name)) +
  geom_point(data = . %>% filter(OR_name != 'other'), aes(alpha = OR_name), size = 2) +
  ggrepel::geom_label_repel(data = . %>% filter(gene_name %in% c('Olfr545')),
                            aes(label = OR_name), force = 0, nudge_y = -1, nudge_x = -0.1, min.segment.length = 0) +
  ggrepel::geom_label_repel(data = . %>% filter(gene_name %in% c('Olfr151')),
                            aes(label = OR_name), force = 0, nudge_y = 4, nudge_x = 0.05, min.segment.length = 0) +
  ggrepel::geom_label_repel(data = . %>% filter(gene_name %in% c('Olfr17')),
                            aes(label = OR_name), force = 0, nudge_y = -3, nudge_x = -0.025, min.segment.length = 0) +
  ggrepel::geom_label_repel(data = . %>% filter(gene_name %in% c('Olfr1507','Olfr16')),
                            aes(label = OR_name), force = 2, nudge_y = 3, nudge_x = -0.05, min.segment.length = 0) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none',
        panel.grid.minor.x = element_blank()) +
  ylab('Stress Score') +
  scale_color_manual(values = new_color_scale_2_17) + 
  scale_alpha_manual(values = new_alpha_scale_2_17) +
  scale_x_continuous(breaks = c(0)) +
  coord_cartesian(ylim = c(-1,1)*max(abs(res[str_detect(res$gene_name,'Olfr'),]$log2FoldChange), na.rm = T))

#Colored for Bright/Dim
set.seed(1)
ggplot(data = res %>%
         filter(str_detect(gene_name, 'Olfr') & !is.na(log2FoldChange)) %>%
         mutate('kde_val' = map_dbl(log2FoldChange, .f = function(x) {
           kde_est <- kde$y[which.min(abs(kde$x - x))]
           0 + runif(1,-kde_est, kde_est)
         })) %>%
         mutate('identity' = ifelse(log2FoldChange >0, 'iRFP_Bright','iRFP_Dim')), 
       aes(x = kde_val, y = log2FoldChange, color = identity)) +
  geom_polygon(data = data.frame('log2FoldChange' = kde$x, 'kde_val' = kde$y) %>%
                 bind_rows(mutate(., 'kde_val' = -kde_val) %>% arrange(-log2FoldChange)), 
               fill = NA, color = 'black', alpha = 1) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_label(data = . %>% 
               group_by(identity) %>%
               summarize() %>%
               mutate('vjust' = ifelse(identity == 'iRFP_Bright', 1, 0),
                      'log2FoldChange' = (vjust-0.5)*Inf,
                      'lab' = fct_recode(identity, 'Higher Stress' = 'iRFP_Bright',
                                         'Lower Stress' = 'iRFP_Dim')),
             aes(label = lab, vjust = vjust), x = -Inf, hjust = 0) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none',
        panel.grid.minor.x = element_blank()) +
  ylab('Stress Score') +
  scale_x_continuous(breaks = c(0)) +
  coord_cartesian(ylim = c(-1,1)*max(abs(res[str_detect(res$gene_name,'Olfr'),]$log2FoldChange), na.rm = T))

#PCA Plot
plotPCA(vst(dds), intgroup = c('Population','Intensity'), returnData = T) %>%
{
  ggplot(data = as.data.frame(.),
         aes(x = PC1, y = PC2, color = Intensity, shape = Population)) +
    geom_point() +
    theme_bw() +
    xlab(paste0('PC1- ', round(attr(.,'percentVar')[1]*100,0),'% Variance')) +
    ylab(paste0('PC2- ', round(attr(.,'percentVar')[2]*100,0),'% Variance')) +
    theme(legend.position = 'bottom', legend.title = element_blank(),
          legend.box = 'vertical', legend.margin=margin()) +
    guides('color' = guide_legend(nrow = 1), 'shape' = guide_legend(nrow = 1)) 
}

#MA Plot of All Genes
ggplot(res %>%
         filter(!is.na(significance)) %>% 
         mutate('identity' = ifelse(str_detect(gene_name, 'Olfr'), 'Olfactory Receptor', 'Other')),
       aes(x = baseMean, y = log2FoldChange, color = significance)) +
  geom_point(alpha = 0.1) +
  facet_wrap(facets = vars(identity)) +
  theme_bw() +
  scale_x_log10() +
  #  geom_text(data = . %>% group_by(identity) %>%
  #              arrange(padj, .by_group = T) %>% dplyr::slice(1:50),
  #            aes(label = gene_name), size = 3, color = 'black') +
  geom_text(data = . %>% filter(!is.na(significance)) %>% 
              group_by(identity, significance) %>% summarize(n=n()) %>%
              mutate('label' = paste0('n=',n),
                     'vjust' = seq(1, length(n))),
            aes(label = label, vjust = vjust), x = Inf, y = Inf, hjust = 1, show.legend = F) +
  geom_text(data = . %>% 
              group_by(identity, significance) %>% summarize(n=n()) %>%
              dcast(significance ~ identity) %>% column_to_rownames(var = 'significance') %>% as.matrix() %>% 
              {tibble('identity' = 'Olfactory Receptor', 'pval' = chisq.test(.)$p.value)},
            aes(label = paste0('p=', formatC(pval, format = 'e', digits = 2))), x = -Inf, y = Inf, hjust = 0, vjust = 1, color = 'black') +
  scale_color_manual(values = c('Significant' = 'red', 'Not Significant' = 'black')) +
  theme(legend.position = 'bottom', legend.title = element_blank()) + 
  guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  xlab('Mean of Normalized Counts') + ylab('Log2 Fold Change\niRFP Bright/Dim') 

#Violin Plot ORs by Zone
ggplot(data = res %>% 
         left_join(zonal_annotation, by = 'gene_name') %>%
         filter(!is.na(zone)) %>%
         mutate('class' = ifelse(zone == 'fishOR', 'Class I','Class II'),
                'zone' = ifelse(zone == 'fishOR','1',as.character(zone))),
       aes(x = zone, y = log2FoldChange, fill = zone)) +
  geom_violin() +
  geom_boxplot(fill = 'white', width = 0.2) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) +
  geom_text(data = . %>% 
  {lm(data = ., formula = log2FoldChange ~ interaction(class,zone))} %>%
  {tibble('pval' = anova(.)$`Pr(>F)`[1],
          'class' = 'Class II')},
  aes(label = paste0('ANOVA p=', formatC(pval, format = 'e', digits = 2))), x = Inf, y = Inf, hjust = 1, vjust = 1, color = 'black', 
  inherit.aes = F) +
  #  geom_label(data = . %>% right_join(data.frame('gene_name' = c('Olfr151','Olfr16','Olfr17', 'Olfr1507'),
  #                                                'OR_name' = c('M71','Mor23','P2','Mor28')))  %>%
  #               mutate('vjust' = ifelse(gene_name %in% c('Olfr17','Olfr16'), 2, -1)), 
  #             aes(vjust = vjust, label = OR_name), fill = 'white') +
  #  geom_point(data = . %>% right_join(data.frame('gene_name' = c('Olfr151','Olfr16','Olfr17', 'Olfr1507'),
  #                                                'OR_name' = c('M71','Mor23','P2','Mor28'))),
  #             color = 'black', fill = 'red', shape = 23, size = 4) +
  facet_grid(cols = vars(class), scales = 'free_x', space = 'free_x') +
  ylab('Stress Score') + guides(fill = guide_legend(nrow = 1))

#ECDF Stress Scores Atf5-Dept vs. Indept
ggplot(data = atf5ko_vs_wt %>%
         dplyr::select(gene_name = 'Gene.Name', KO = 'AtfeKO-OmpGFP+1', WT ='WT-OMPGFP+1') %>%
         mutate('log2FC_KO_vs_WT' = log2(KO/WT),
                'identity' = ifelse(log2FC_KO_vs_WT > 0, 'Atf5-Independent','Atf5-Dependent')) %>%
         right_join(res, by = 'gene_name') %>%
         filter(str_detect(gene_name, 'Olfr') & !is.na(identity)),
       aes(x = log2FoldChange, color = identity)) +
  stat_ecdf() +
  geom_text(data = . %>% group_by(identity) %>% summarize(n=n()),
            aes(label = paste0('n=', n), vjust = as.numeric(as.factor(identity))), x = -Inf, y = Inf, hjust = 0, show.legend = F) +
  theme_bw() +
  # ggtitle('Atf5 Dependence in Sorted Omp-GFP Cells') +
  xlab('Stress Score') + ylab('Cumulative Density') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2))

#Lack of Zonal Effects Atf5 KO vs. WT and Rtp KO vs. WT
ggplot(data = atf5ko_vs_wt %>%
         dplyr::select(gene_name = 'Gene.Name', KO = 'AtfeKO-OmpGFP+1', WT ='WT-OMPGFP+1') %>%
         mutate('log2FC_KO_vs_WT' = log2(KO/WT),
                'identity' = ifelse(log2FC_KO_vs_WT > 0, 'Atf5-Independent','Atf5-Dependent')) %>%
         left_join(zonal_annotation, by = 'gene_name') %>%
         filter(!is.na(zone)) %>%
         mutate('facet' = ifelse(zone == 'fishOR', 'Class I','Class II'),
                'zone' = ifelse(zone == 'fishOR','1',as.character(zone))),
       aes(x = zone, y = log2FC_KO_vs_WT)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_violin(aes(fill = zone)) +
  geom_boxplot(fill = 'white', width = 0.2, outlier.shape = NA) +
  theme_bw() +
  facet_grid(cols = vars(facet), scales = 'free_x', space = 'free_x') +
  ggtitle('Atf5 Dependence') + ylab('Log2FC Atf5 KO vs. WT') +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), legend.title = element_blank(), 
        legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 1))

#ECDF Stress Score ~ Rtp Dependence
ggplot(data = res %>% 
         mutate('rtp_dependance' = map_chr(gene_name, .f = function(x) {
           if (x %in% matsunami_ORs$oOR) {'Rtp-Independent'}
           else if (x %in% matsunami_ORs$uOR) {'Rtp-Dependent'}
           else {'Other'}
         })) %>% 
         filter(str_detect(gene_name, 'Olfr') & !is.na(rtp_dependance) & (rtp_dependance != 'Other')),
       aes(x = log2FoldChange, color = rtp_dependance)) +
  stat_ecdf() +
  geom_text(data = . %>% group_by(rtp_dependance) %>% summarize(n=n()),
            aes(label = paste0('n=', n), vjust = as.numeric(as.factor(rtp_dependance))), x = -Inf, y = Inf, hjust = 0, show.legend = F) +
  theme_bw() +
  # ggtitle('Rtp1/2 Dependence in Sorted Omp-GFP Cells') +
  xlab('Stress Score') + ylab('Cumulative Density') +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  guides(color = guide_legend(nrow = 2))

#Show that Rtp1/2 Dependence has NO zonal confound
ggplot(data = matsunami_de %>%
         filter(str_detect(gene_name,'Olfr')) %>%
         mutate('rtp_dependance' = map_chr(gene_name, .f = function(x) {
           if (x %in% matsunami_ORs$oOR) {'Rtp-Independant'}
           else if (x %in% matsunami_ORs$uOR) {'Rtp-Dependant'}
           else {'Other'}
         })) %>%
         left_join(zonal_annotation, by = 'gene_name') %>%
         filter(!is.na(zone)) %>%
         mutate('facet' = ifelse(zone == 'fishOR', 'Class I','Class II'),
                'zone' = ifelse(zone == 'fishOR','1', as.character(zone))),
       aes(x = zone, y = log2FoldChange)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_violin(aes(fill = zone)) +
  geom_boxplot(fill = 'white', width = 0.2, outlier.shape = NA) +
  facet_grid(cols = vars(facet), scales = 'free_x', space = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), legend.title = element_blank(), 
        legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 1)) +
  ggtitle('Rtp1/2 Dependence') + ylab('Log2FC Rtp KO vs. WT') 

#Stress Scores in mOSNs ONLY
dds_omponly <- DESeqDataSetFromTximport(txi = tximport(file_paths[coldata$Population == 'iRFPGFP'], 
                                                       type = 'salmon', tx2gene = tx2gene), 
                                        colData = coldata[coldata$Population == 'iRFPGFP',],
                                        design = ~Intensity) %>%
  DESeq()

res_omponly <- results(dds_omponly, contrast = c('Intensity','High','Low')) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  mutate('significance' = ifelse(padj < 0.05, 'Significant','Not Significant'))

#Simpler version of this plot for Supplemental Figure
ggplot(data = res %>%
         filter(str_detect(gene_name,'Olfr')) %>%
         left_join(res_omponly, by = c('gene_name'),
                   suffix = c('_all','_omp')) %>%
         filter(!(is.na(log2FoldChange_all) | is.na(log2FoldChange_omp))),
       aes(x = log2FoldChange_all, y = log2FoldChange_omp)) +
  geom_point(alpha = 0.1) +
  theme_bw() +
  geom_text(data = . %>% 
              summarize(rho=cor.test(log2FoldChange_all, log2FoldChange_omp, method = 'spearman')$estimate, .groups = 'drop'),
            aes(label = paste0('rho=', round(rho,2))), x = -Inf, y = Inf, vjust = 1, hjust = 0) +
  stat_function(fun = function(x) {x}, linetype = 'dashed', color = 'black') +
  geom_smooth(method = 'lm', se = T, alpha = 0.2, color = 'blue') +
  xlab('Stress Score- iOSNs and mOSNs') + ylab('Stress Score- mOSNs only')