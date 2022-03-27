#RNA-seq on developing OSNs

library(tidyverse)
library(reshape2)
library(DESeq2)
library(patchwork)

#Load Files
salmon_ls <- list.files('./',
                        pattern = '.*quant.sf', recursive = T, full.names = T) %>%
  .[str_detect(.,'INP|iOSN|mOSN')] %>% 
  set_names(map_chr(.,.f = function(x) {x %>% str_split('/') %>% .[[1]] %>% .[length(.)-1] %>% 
      str_extract('(?<=SalmonOut_).*')}))

tx2gene <- read.table('tx_to_gene_name.tsv', sep='\t', 
                      header=F, stringsAsFactors = F, col.names = c('tx','gene'))

#Import for DESeq
imported_data <- tximport::tximport(salmon_ls, type = 'salmon', tx2gene = tx2gene) 

#ColData
salmon_coldata <- data.frame('lib_id' = names(salmon_ls)) %>%
  mutate('pop' = str_extract(lib_id, 'INP|iOSN|mOSN'),
         'pop' = fct_relevel(pop, 'INP','iOSN')) %>%
  column_to_rownames(var = 'lib_id')

#DESeq Dataset
dds_salmon <- DESeqDataSetFromTximport(txi = imported_data, 
                                       colData = salmon_coldata,
                                       design = ~pop) %>%
  DESeq(test = 'LRT',reduced = ~1)

pca_data <- plotPCA(vst(dds_salmon), intgroup = c('pop'), returnData = T)
ggplot(data = pca_data,
       aes(x = PC1, y = PC2, color = pop)) +
  geom_point() +
  theme_bw() +
  xlab(paste0('PC1-', round(attr(pca_data,'percentVar')[1]*100,0),'% Variance')) +
  ylab(paste0('PC2-', round(attr(pca_data,'percentVar')[2]*100,0),'% Variance')) +
  theme(legend.title = element_blank())

#Extract the LRT Results from FACS sorted cells
res_differentiation_lrt <- map_dfr(resultsNames(dds_salmon)[-1], .f = function(x) {
  results(dds_salmon, contrast = list(x)) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene_name') %>%
    mutate('contrast' = rep(x, nrow(.)))
}) %>%
  dcast(gene_name + baseMean + stat + pvalue + padj ~ contrast, value.var = 'log2FoldChange') 

#Prepare a heatmap for significant genes
n_clusters <- 6
vst_mat_all <- res_differentiation_lrt %>%
  filter(padj < 0.05) %>% 
  pull(gene_name) %>%
  {assays(vst(dds_salmon))[[1]][.,]}

row_dendrogram <- hclust(dist(t(scale(t(vst_mat_all))))) #dendrogram
group_membership <- cutree(row_dendrogram, k = n_clusters) %>% as.character() #arbitrary groups, will reorder below
plot_order <- row_dendrogram$order #order of rows in vst_mat_all/group_membership to plot in dendrogram order

group_order <- seq(1,n_clusters) %>% as.character() #order we want for the facetted plots (REVERSE of dendrogram...)

#Rename the arbitrary groups so that they go 1..6 in final figure (6...1 for dendrogram)
lookup_table <- 
  data.frame('input_group_order' = rev(unique(group_membership[plot_order])) %>% as.character(), #names of group, in REVERSE order encountered when plotting dendrogram
             'out_order' = seq(1, length(unique(group_membership)))) #order we want

group_membership <- 
  data.frame('input_group_order' = group_membership) %>%
  left_join(lookup_table, by = 'input_group_order') %>%
  pull(out_order)

col_order <- colData(dds_salmon) %>% 
  as.data.frame() %>%
  rownames_to_column(var = 'lib_name') %>% 
  arrange(as.numeric(pop), lib_name) %>%
  pull(lib_name)

#Heatmap, clustered 
heatmap_lrt <- 
  ggplot(data = as.data.frame(t(scale(t(vst_mat_all)))) %>%
           rownames_to_column(var = 'gene_name') %>%
           mutate('gene_name' = factor(gene_name, levels = gene_name[plot_order])) %>%
           melt(id.vars = 'gene_name') %>%
           mutate('variable' = factor(variable, levels = col_order)),
         aes(x = variable, y = gene_name, fill = value)) +
  geom_tile() +
  theme_void() + theme(legend.position = 'none') +
  scale_fill_gradient2(low = 'red',mid = 'black',high = 'green') +
  scale_x_discrete(labels = function(x) {
    lib_type <- str_extract(x,'Atf5|Ngn|Omp')
    rep_type <- str_extract(x,'[0-9]$')
    paste0(lib_type, '_', rep_type)
  }) +
  xlab('') + ylab('') 

#Clusters from Cut Dendrogram
annotation_lrt <- 
  ggplot(data = 
           data.frame('gene_name' = rownames(vst_mat_all),
                      'cluster' = group_membership) %>%
           mutate('cluster' = factor(group_membership, levels = group_order),
                  'gene_name' = factor(gene_name, levels = gene_name[plot_order])),
         aes(x = 1, y = gene_name, fill = cluster)) +
  geom_tile() +
  theme_void() +
  theme(legend.position = 'none') +
  geom_rect(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, color = 'black', fill = NA)

#Summary of Clusters
summary_plot <- 
  ggplot(data = as.data.frame(t(scale(t(vst_mat_all)))) %>%
           mutate('group_membership' = group_membership) %>% 
           melt(id.vars = c('group_membership')) %>%
           left_join(rownames_to_column(as.data.frame(colData(dds_salmon)), var = 'variable')) %>%
           group_by(group_membership, pop) %>%
           summarize(median = median(value),
                     'q1' = quantile(value, 0.25),
                     'q3' = quantile(value, 0.75)) %>%
           mutate('group_membership' = factor(group_membership, group_order),
                  'cell_type' = fct_relevel(pop,'INP','iOSN')), # for now
         aes(x = cell_type, y = median)) +
  geom_point() + 
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0) +
  geom_line(aes(group = group_membership)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_grid(rows = vars(group_membership)) +
  geom_rect(data = . %>% group_by(group_membership) %>% summarize(),
            aes(fill = group_membership), 
            color = 'black', xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3, inherit.aes = F,
            show.legend = F) +
  theme_bw() +
  theme(legend.position = 'bottom', axis.title.x = element_blank()) +
  ylab('Scaled Expression')


#Add the GO Analysis. No ORs
go_out <- 
{{map(group_order, .f = function(x) {
  rownames(vst_mat_all)[group_membership == x] %>%
    .[!str_detect(.,'Olfr')]
}) %>%
    gprofiler2::gost(organism = 'mmusculus', multi_query = T, 
                     custom_bg = res_differentiation_lrt$gene_name %>%
                       .[!str_detect(.,'Olfr')],
                     sources = 'GO') %>%
    riboHelpers::plot_gProfiler(clusters_vector = group_order, 
                                regex_highlight = 'ER|endoplasmic|unfolded|stress|topolo',
                                max_categories = 30)} %>%
    set_names(group_order) %>% 
    imap(.f = function(x,n) {x + ggtitle(n) + theme(plot.title = element_text(hjust = 1))}) %>%
    purrr::reduce(.f = `+`, .dir = 'forward')} +
  plot_layout(ncol = 1)

#Add just 2 GO analysis for clusters 5 and 6, No OR's
clusters_to_highlight <- c('5','6')

gprofiler_out_to_plot <- 
  map(clusters_to_highlight, .f = function(x) {
    rownames(vst_mat_all)[group_membership == x] %>%
      .[!str_detect(.,'Olfr')]
  }) %>%
  gprofiler2::gost(organism = 'mmusculus', multi_query = T, 
                   custom_bg = res_differentiation_lrt$gene_name %>%
                     .[!str_detect(.,'Olfr')],
                   sources = 'GO',
                   correction_method = 'fdr') %>%
                   {.$result} %>%
  dplyr::select(-parents) %>%
  unnest() %>% 
  group_by(term_id, term_size, source, term_name) %>% 
  mutate(cluster = clusters_to_highlight) %>% ungroup() %>%
  filter(p_values < 0.05) %>%
  group_by(cluster) %>%
  arrange(p_values,.by_group = T) %>%
  dplyr::slice(1:15) %>%
  arrange(desc(p_values), .by_group = T) %>%
  mutate('term_name_to_plot' = paste0(cluster,'.',term_name), #avoid overlap of term appearing for >1 cluster
         'term_name_to_plot' = factor(term_name_to_plot, levels = term_name_to_plot))

go_out_5_6 <- 
  ggplot(data = gprofiler_out_to_plot,
         aes(x = term_name_to_plot, y = -log10(p_values))) +
  geom_rect(data = . %>% group_by(cluster) %>%
              summarize(),
            aes(fill = cluster), inherit.aes = F, 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
  geom_col(fill = 'green', color = 'black') +
  theme_bw() +
  coord_flip() +
  facet_grid(rows = vars(cluster), scales = 'free_y', space = 'free_y') +
  scale_x_discrete(labels = function(x) {
    str_wrap(str_extract(x,'(?<=\\.).*'),50)}) +
  ylab('-log10(FDR)') + xlab('') + theme(legend.position = 'none') +
  scale_fill_manual(values = scales::hue_pal()(6) %>% set_names(c(1:6)), drop = F)

#Plot Everything
(heatmap_lrt + annotation_lrt + summary_plot + go_out_5_6 + 
    plot_layout(widths = c(0.4,0.05,0.2,0.3)))

#Look at counts for all Olfrs
vst_all_olfr <- rownames(dds_salmon) %>%
  .[str_detect(.,'Olfr')] %>%
  {assays(vst(dds_salmon))[[1]][.,]}

#Generate the Z-Scored Plot
ggplot(data = t(scale(t(vst_all_olfr))) %>%
         as.data.frame() %>%
         rownames_to_column(var = 'gene_name') %>%
         melt(id.vars = c('gene_name')) %>%
         left_join(rownames_to_column(as.data.frame(colData(dds_salmon)), var = 'variable'), by = 'variable') %>%
         group_by(gene_name, pop) %>% 
         summarize('value' = mean(value)) %>%
         ungroup() %>%
         mutate('pop' = fct_relevel(pop, 'INP','iOSN')),
       aes(x = pop, y = value)) +
  geom_line(aes(group = gene_name), alpha = 0.01) +
  geom_point(data = . %>%
               group_by(pop) %>%
               summarize(value = median(value, na.rm = T)),
             color = 'red') +
  geom_line(data = . %>%
              group_by(pop) %>%
              summarize(value = median(value, na.rm = T)),
            color = 'red', group = 1) +
  geom_errorbar(data = . %>%
                  group_by(pop) %>%
                  summarize('q1' = quantile(value, 0.25, na.rm = T),
                            'q3' = quantile(value, 0.75, na.rm = T),
                            value = median(value, na.rm = T)),
                aes(ymax = q3, ymin = q1), color = 'red', width = 0) +
  theme_bw() +
  ylab('Scaled OR Expression') + xlab('')

