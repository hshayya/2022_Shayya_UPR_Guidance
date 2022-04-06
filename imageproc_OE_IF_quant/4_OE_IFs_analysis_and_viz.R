library(tidyverse)
library(reshape2)
library(dbscan)
setwd('/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Experiments/Analysis_Glomeruli/IF_Quants/OE_IF_Quants/')

#Semi-Quantiative, local, analysis of Mor28 IF signal intensities in GFP and tdtom cells in M28 Ddit3 WT vs. cKO mice
#Cell locations/mean M28 IF signal read into R -> hdbscan clustering to create local environments -> measure avg M28 in tdtom - GFP cells within each local environment -> plot

#Read Cell Measurements
cell_measurements <- read_tsv('cell_measurements.tsv') %>%
  mutate('geno' = str_extract(cell, 'WT|Exp'), #WT= WT, Exp = cKO
         'section' = str_extract(cell,'.*(?=_tdtom|_gfp)'),
         'cell_type' = str_extract(cell, 'tdtom|gfp')) %>%
  mutate('cell_type' = fct_relevel(cell_type, 'gfp'),
         'geno' = fct_relevel(geno,'WT'))

#Read Cell Positions
cell_positions <- read_tsv('cell_annotations.tsv')

#For Plotting
gfp_tdtom_mapping <- c('gfp' = 'green','tdtom' = 'red')

#Cutoffs for Cell Size
ggplot(data = cell_measurements,
       aes(x = area, group = cell_type, color = cell_type)) +
  geom_line(stat = 'density') +
  facet_grid(rows = vars(geno)) +
  theme_bw() +
  geom_vline(data = data.frame('area' = c(30,100)), 
             aes(xintercept = area), linetype = 'dashed') +
  scale_color_manual(values = gfp_tdtom_mapping)

#Clustering Approach- Find hdbscan > kmeans. 
set.seed(2)
clustered_cells_hdb <- cell_measurements %>% 
  left_join(cell_positions, by = 'cell') %>%
  group_by(section) %>%
  mutate('cluster' = hdbscan(as.matrix(data.frame('x_um' = x_um, 'y_um' = y_um)), 
                             minPts = 6)$cluster) %>%
  ungroup() 

#Verify that local environments look good
ggplot(data = clustered_cells_hdb %>%
         mutate('cluster' = ifelse(cluster == 0, NA,cluster),
                'cluster' = factor(cluster)),
       aes(x = x_um, y = -y_um, color = cluster, shape = cell_type)) +
  geom_point() +
  facet_wrap(facets = vars(section), nrow = 2) +
  theme_bw()

#Mean log10-M28 signal GFP vs. tdtom cells ~Genotype (each point = 1 local environment)
ggplot(data = clustered_cells_hdb %>%
         filter(cluster != 0) %>% 
         mutate('log10_mean' = log10(mean)) %>% #should always look at this in log space
         group_by(section, cell_type, cluster, geno) %>%
         summarize('log10_mean' = mean(log10_mean),
                   'n' = n()) %>%
         group_by(section, geno, cluster) %>%
         summarize('log10_change_mean' = log10_mean[cell_type == 'tdtom'] - log10_mean[cell_type == 'gfp'],
                   'n_total' = sum(n)) %>%
         ungroup() %>%
         mutate('geno' = fct_recode(geno,'cKO' = 'Exp')),
       aes(x = geno, y = log10_change_mean, color = geno)) +
  geom_point(position = position_jitter(width = 0.1)) +
  geom_boxplot(width = 0.2, fill = NA) +
  theme_bw() +
  ylab('Per Cluster M28 Signal in\nTdt - GFP Cells') +
  theme(legend.position = 'none') + xlab('') +
  scale_color_manual(values = scales::hue_pal()(3)[c(1,3)] %>% set_names(c('WT','cKO'))) #match WT/Het/cKO color in other figures

#P value
clustered_cells_hdb %>%
  filter(cluster != 0) %>% 
  mutate('log10_mean' = log10(mean)) %>% #should always look at this in log space
  group_by(section, cell_type, cluster, geno) %>%
  summarize('log10_mean' = mean(log10_mean),
            'n' = n()) %>%
  group_by(section, geno, cluster) %>%
  summarize('log10_change_mean' = log10_mean[cell_type == 'tdtom'] - log10_mean[cell_type == 'gfp'],
            'n_total' = sum(n)) %>%
  ungroup() %>%
  mutate('geno' = fct_recode(geno,'cKO' = 'Exp')) %>%
  summarize('p' = wilcox.test(log10_change_mean[geno == 'cKO'],
                              log10_change_mean[geno == 'WT'])$p.value)