#Calculate Stress Scores in M71 -> P2 Swap OSNs

library(tidyverse)
library(reshape2)
library(DESeq2)
library(tximport)

#Load Files
file_paths <- list.files(path = './',
                         pattern = '.*quant.sf', recursive = T) %>%
  .[str_detect(.,'Ctrl|Swap')] %>%
  set_names(str_extract(.,'.*(?=_quant.sf$)'))

tx2gene <- read.table('tx_to_gene_name_gfp_tdtom_lacz.tsv', sep='\t', 
                      header=F, stringsAsFactors = F, col.names = c('tx','gene'))

zonal_annotation <- read.table('ORs-by-zone.txt', sep = '\t', header = F,
                               col.names = c('gene_name', 'zone'))

#Prep for DESeq2
imported_data <- tximport(file_paths, type = 'salmon', tx2gene = tx2gene)

coldata <- data.frame('lib_id' = names(file_paths)) %>%
  mutate('Group' = str_extract(lib_id, 'Ctrl|Swap'),
         'Animal' = str_extract(lib_id, '(?<=Ctrl|Swap)[0-9]+'),
         'Population' = str_extract(lib_id,'iRFP(GFP)?'),
         'Intensity' = str_extract(lib_id,'Low|High'),
         'Group' = fct_relevel(Group, 'Ctrl','Swap'),
         'Population' = fct_relevel(Population,'iRFP','iRFPGFP'),
         'Intensity' = fct_relevel(Intensity, 'Low','High')) %>%
  column_to_rownames(var = 'lib_id')

#Check Full Dataset first for clustering and the like
dds_full <- DESeqDataSetFromTximport(txi = imported_data, 
                                     colData = coldata,
                                     design = ~Population+Intensity) %>%
  DESeq()

#PCA for Swap Experiment
pca_obj_swap_final <- plotPCA(vst(dds_full), 
                              intgroup = c('Population','Intensity', 'Group'), returnData = T)
ggplot(data = pca_obj_swap_final,
       aes(x = PC1, y = PC2, color = interaction(Population, Intensity, lex.order = T),
           shape = Group)) +
  geom_point() +
  theme_bw() +
  xlab(paste0('PC1-', round(attr(pca_obj_swap_final, 'percentVar')[1]*100), '% Variance')) +
  ylab(paste0('PC2-', round(attr(pca_obj_swap_final, 'percentVar')[2]*100), '% Variance')) +
  theme(legend.position = 'bottom', legend.box = 'horizontal', legend.title = element_blank()) +
  guides('color' = guide_legend(nrow = 2, order = 1),
         'shape' = guide_legend(nrow = 2, order = 2))

#Compute Final Stress Scores
stress_scores_per_animal_final <- counts(dds_full, normalized = T) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  melt(id.vars = c('gene_name'),
       variable.name = 'lib_id') %>%
  left_join(rownames_to_column(as.data.frame(colData(dds_full)), var = 'lib_id'),
            by = 'lib_id') %>%
  filter(str_detect(gene_name,'Olfr|lacZ')) %>%
  group_by(gene_name, Group, Animal, Population) %>%
  summarize('log2' = log2(value[Intensity == 'High']/value[Intensity == 'Low'])) %>%
  summarize('log2FoldChange' = mean(log2, na.rm = T)) %>%
  ungroup()

#write_tsv(stress_scores_per_animal_final,'Swap_stressscores_peranimal.tsv')

#Point Plot for M71/P2/Swap, For Paper
new_color_scale_2_17 <- scales::hue_pal()(7)[c(1,3,5,6,7)] %>% set_names(c('M71','Class I','Mor23','P2','Mor28')) %>%
  c('other' = 'black')
ggplot(stress_scores_per_animal_final %>%
         filter(gene_name %in% c('lacZ','Olfr151','Olfr17')) %>%
         mutate('gene_name' = fct_recode(gene_name,'M71'='Olfr151','P2' = 'Olfr17','M71->P2' = 'lacZ'),
                'gene_name' = fct_relevel(gene_name, 'M71','P2','M71->P2')) %>%
         filter(is.finite(log2FoldChange) & 
                  ((Group == 'Swap' & gene_name == 'M71->P2') | 
                     (Group == 'Ctrl')
                  )),
       aes(x = gene_name, y = log2FoldChange, color = gene_name)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c(new_color_scale_2_17, 'M71->P2' = as.vector(new_color_scale_2_17)[2]), 
                     na.value = 'black') +
  xlab('') + ylab('Stress Score') + 
  theme(legend.position = 'none')

#P value, Final, For Paper
stress_scores_per_animal_final %>%
  filter(gene_name %in% c('lacZ','Olfr151','Olfr17')) %>%
  mutate('gene_name' = fct_recode(gene_name,'M71'='Olfr151','P2' = 'Olfr17','M71->P2' = 'lacZ'),
         'gene_name' = fct_relevel(gene_name, 'M71','P2','M71->P2')) %>%
  filter(is.finite(log2FoldChange) & 
           ((Group == 'Swap' & gene_name == 'M71->P2') | 
              (Group == 'Ctrl')
           )) %>%
  lm(formula = log2FoldChange ~ gene_name) %>%
  aov() %>%
  TukeyHSD()

#Stress Scores in n=6 WT vs. n=3 Swap Libraries
#Note this is NOT per animal -> FIRST taking averages per OR in each group, THEN log2FCs
#We only needed the per-animal approach so that we could show the M71/P2 and LacZ on same scale w/ some statistical measure between them
ggplot(data = counts(dds_full, normalized = T) %>%
         as.data.frame() %>%
         rownames_to_column(var = 'gene_name') %>%
         melt(id.vars = c('gene_name'),
              variable.name = 'lib_id') %>%
         left_join(rownames_to_column(as.data.frame(colData(dds_full)), var = 'lib_id'),
                   by = 'lib_id') %>%
         filter(str_detect(gene_name,'Olfr')) %>%
         group_by(gene_name, Group, Population, Intensity) %>%
         summarize('mean_cts' = mean(value)) %>%
         summarize('log2' = log2(mean_cts[Intensity == 'High']/mean_cts[Intensity == 'Low'])) %>%
         summarize('log2FoldChange' = mean(log2, na.rm = T)) %>%
         dcast(gene_name ~ Group, value.var = 'log2FoldChange') %>%
         filter(complete.cases(.) & is.finite(Ctrl) & is.finite(Swap)),
       aes(x = Ctrl, y = Swap)) +
  geom_point(alpha = 0.1) +
  theme_bw() +
  stat_function(fun = function(x) {x}, linetype = 'dashed', color = 'black') +
  geom_smooth(method = 'lm', se = T, color = 'blue', alpha = 0.2) +
  geom_text(data = . %>% summarize('est' = cor.test(Ctrl,Swap, method = 'spearman')$estimate),
            aes(label = paste0('rho=', round(est,2)), x = -Inf, y = Inf, hjust = 0, vjust = 1)) +
  xlab('Atf5(rep/+); Omp(iGFP/+)') + ylab('Atf5(rep/+); Omp(iGFP/+);\nP2(M71iLacZ/+)')

#Sina Plot
#WT(n=6), Swap(n=3) -> most genes n=9 except M71&P2(WT only, n=6) and lacZ (Swap only, n=3)
#These values ONLY used for swap experiment. 
#For all other experiments use n=4 defined from the below DESeq2 approach
stress_scores_for_swap_sina <- counts(dds_full, normalized = T) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  melt(id.vars = c('gene_name'),
       variable.name = 'lib_id') %>%
  left_join(rownames_to_column(as.data.frame(colData(dds_full)), var = 'lib_id'),
            by = 'lib_id') %>%
  filter(str_detect(gene_name,'Olfr|lacZ')) %>%
  filter(!(gene_name == 'lacZ' & Group == 'Ctrl') & 
           !(gene_name %in% c('Olfr151','Olfr17') & Group == 'Swap')) %>%
  group_by(gene_name, Population, Intensity) %>%
  summarize('mean_cts' = mean(value)) %>%
  summarize('log2' = log2(mean_cts[Intensity == 'High']/mean_cts[Intensity == 'Low'])) %>%
  summarize('log2FoldChange' = mean(log2, na.rm = T)) %>%
  ungroup()

#write_tsv(stress_scores_for_swap_sina, 'Swap_stressscores.tsv')

#KDE for Swap Sina Plot
kde_swap_sina <- stress_scores_for_swap_sina %>%
  filter(!is.na(log2FoldChange) & is.finite(log2FoldChange))  %>%
  {density(.$log2FoldChange)}

#Sina Plot To find optimal seed, mapped across seeds looking for all kde_val < 0.05
new_alpha_scale_2_17 <- ifelse(names(new_color_scale_2_17) == 'other', 0.1, 1) %>% set_names(names(new_color_scale_2_17))

set.seed(1317)
ggplot(data = stress_scores_for_swap_sina %>%
         filter(!is.na(log2FoldChange) & is.finite(log2FoldChange)) %>%
         mutate('kde_val' = map_dbl(log2FoldChange, .f = function(x) {
           kde_est <- kde_swap_sina$y[which.min(abs(kde_swap_sina$x - x))]
           0 + runif(1,-kde_est, kde_est)
         })) %>%
         left_join(zonal_annotation, by = 'gene_name') %>%
         left_join(data.frame('gene_name' = c('Olfr151','Olfr16','Olfr17', 'Olfr1507','Olfr545', 'lacZ'),
                              'OR_name' = c('M71','Mor23','P2','Mor28','Class I', 'M71->P2'), stringsAsFactors = F)) %>%
         mutate('OR_name' = replace_na(OR_name, 'other'),
                'zone' = replace_na(as.character(zone), 'other'),
                'OR_name' = ifelse(OR_name == 'other', ifelse(zone == 'fishOR','Class I','other'), OR_name),
                'OR_name' = fct_relevel(OR_name,'M71->P2','Class I','M71','Mor23','P2','Mor28','other'),
                'OR_name' = fct_relevel(OR_name, rev)), 
       aes(x = kde_val, y = log2FoldChange, color = OR_name)) +
  geom_polygon(data = data.frame('log2FoldChange' = kde_swap_sina$x, 'kde_val' = kde_swap_sina$y) %>%
                 bind_rows(mutate(., 'kde_val' = -kde_val) %>% arrange(-log2FoldChange)), 
               fill = NA, color = 'black', alpha = 1) +
  geom_point(data = . %>% filter(!(OR_name %in% c('M71','P2','M71->P2'))), alpha = 0.075, color = 'black') +
  geom_point(data = . %>% filter(OR_name %in% c('M71','P2','M71->P2')), alpha = 1, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none',
        panel.grid.minor.x = element_blank()) +
  ylab('Stress Score') +
  scale_color_manual(values = c(new_color_scale_2_17, 'M71->P2' = as.vector(new_color_scale_2_17)[2]), 
                     na.value = 'black') + 
  scale_x_continuous(breaks = c(0)) +
  ggrepel::geom_label_repel(data = . %>% filter(gene_name %in% c('Olfr151')),
                            aes(label = OR_name), force = 0, nudge_y = 1.5, nudge_x = -0.1, min.segment.length = 0) +
  ggrepel::geom_label_repel(data = . %>% filter(gene_name %in% c('lacZ')),
                            aes(label = OR_name), force = 0, nudge_y = 1.5, nudge_x = 0.1, min.segment.length = 0) +
  ggrepel::geom_label_repel(data = . %>% filter(gene_name %in% c('Olfr17')),
                            aes(label = OR_name), force = 0, nudge_y = -2.5, nudge_x = 0.1, min.segment.length = 0) +
  coord_cartesian(ylim = c(-1,1)*max(abs(stress_scores_for_swap_sina[is.finite(stress_scores_for_swap_sina$log2FoldChange),]$log2FoldChange), na.rm = T))