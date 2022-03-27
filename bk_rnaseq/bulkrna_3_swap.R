#Calculate Stress Scores in M71 -> P2 Swap OSNs

library(tidyverse)
library(reshape2)
library(DESeq2)
library(tximport)

#Load Files
file_paths <- list.files(path = './',
                         pattern = '.*quant.sf', recursive = T, full.names = TRUE) %>%
  .[str_detect(.,'Ctrl|Swap')] %>%
  set_names(str_extract(.,'(?<=SalmonOut_).*(?=/quant.sf$)'))

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
         'Intensity' = str_extract(lib_id,'Low|High')) %>%
  column_to_rownames(var = 'lib_id')

#Check Full Dataset first for clustering and the like
dds_full <- DESeqDataSetFromTximport(txi = imported_data, 
                                     colData = coldata,
                                     design = ~Group) %>%
  DESeq()

#Compute the Stress scores on a per-lib basis
stress_scores_per_animal <- counts(dds_full, normalized = T) %>%
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

#write_tsv(stress_scores_per_animal,'Swap_stressscores_peranimal.tsv')

#Plot per-animal stress scores for M71, P2 and Swap OSNs
color_scale <- c('Class I','M71','Mor23','P2','Mor28','other')
color_scale <- c(scales::hue_pal()(length(color_scale)-1), 'black') %>% set_names(color_scale)

ggplot(stress_scores_per_animal %>%
         filter(gene_name %in% c('lacZ','Olfr151','Olfr17')) %>%
         mutate('gene_name' = fct_recode(gene_name,'M71'='Olfr151','P2' = 'Olfr17','M71->P2' = 'lacZ'),
                'gene_name' = fct_relevel(gene_name, 'M71','P2','M71->P2')) %>%
         filter(is.finite(log2FoldChange) & 
                  ((Group == 'Swap' & gene_name == 'M71->P2') | 
                     (Group == 'Ctrl') #only consider M71/P2 from control mice (avoid confounding from swap transcript)
                  )),
       aes(x = gene_name, y = log2FoldChange, color = gene_name)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = set_names(color_scale, c('Class I','M71','M71->P2','P2','Mor28','other'))) +
  xlab('') + ylab('Stress Score') + 
  theme(legend.position = 'none')

#Calculate p value for Swap vs. P2 using per-animal scores and Tukey post-test
stress_scores_per_animal %>%
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

#Consistency of stress scores in Atf5(rep/+); Omp(iGFP/+) vs Swap Mice
#Calculate averages within group (Ctrl/Swap, Pop, High/Low), then take average LFCs instead of computing per-animal (lots of dropout)
stress_scores_by_hand <- counts(dds_full, normalized = T) %>%
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
  summarize('log2FoldChange' = mean(log2, na.rm = T))

ggplot(data = stress_scores_by_hand %>%
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

#Sina Plot of Stress Scores for swap Exp 
#This is n=6 for most genes except M71 & P2 (n=4 Ctrls only) and lacZ (n=2 Swap only))
#This is ONLY used for swap experiment. For all other experiments use n=4 defined from the below DESeq2 approach
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

#Scales
#To keep consistency with prior colors from FACS experiments.
new_color_scale_2_17 <- scales::hue_pal()(7)[c(1,3,5,6,7)] %>% set_names(c('M71','Class I','Mor23','P2','Mor28')) %>%
  c('other' = 'black')

new_alpha_scale_2_17 <- ifelse(names(new_color_scale_2_17) == 'other', 0.1, 1) %>% set_names(names(new_color_scale_2_17))

#KDE for Swap Sina Plot
kde_swap_sina <- stress_scores_for_swap_sina %>%
  filter(!is.na(log2FoldChange) & is.finite(log2FoldChange))  %>%
  {density(.$log2FoldChange)}

#Generate the Swap Sina Plot
set.seed(547)
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
