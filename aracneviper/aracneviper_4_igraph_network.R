#Represent Ddit3 Regulon using iGraph

library(igraph)
library(tidyverse)
library(ggraph)
library(reshape2)

#Load Regulons
osn_regulon <- viper::aracne2regulon(afile = 'aracne_network_out_noheader.txt',
                                     eset = '272samples_OE_expression_matrix_forArachne.tsv',
                                     format = '3col')

axon_guidance_genes <- read_tsv('guidance_network.tsv') #see scrna code
scrna_de <- read_tsv('DE_scRNAseq_Bright_Dim_withinzone.tsv') #see scrna code

ddit3_net <- osn_regulon$Ddit3 %>%
{enframe(.[[1]], name = 'Target',value = 'tfmode')} %>%
  mutate('likelihood' = osn_regulon$Ddit3$likelihood,
         'Regulator' = rep('Ddit3',nrow(.)),
         'mode_of_reg' = ifelse(tfmode>0,'Positive','Negative'),
         'guidance_ident' = ifelse(Target %in% axon_guidance_genes$gene_name | 
                                     str_detect(Target,'Ddit3|Rtp|Unc|Cng|Cntn|Eph|Pcdh|Sema|Ten'),
                                   'Guidance Mol','Other')) %>%
  dplyr::select(Regulator, Target, mode_of_reg, likelihood, guidance_ident)

vert_df <- unlist(ddit3_net[,c('Regulator','Target')]) %>%
  unique() %>%
  {data.frame(name = .)} %>%
  left_join(scrna_de, by = c('name' = 'gene_name')) %>%
  mutate('identity' = ifelse(log2FoldChange > 0,'High Stress','Low Stress'),
         'guidance_ident' = ifelse(name %in% axon_guidance_genes$gene_name | 
                                     str_detect(name,'Ddit3|Rtp|Unc|Cng|Cntn|Eph|Pcdh|Sema|Ten'),
                                   'Guidance Mol','Other')) %>%
  dplyr::select(name,identity, guidance_ident)

ddit3_igraph <- graph_from_data_frame(d = ddit3_net,
                                      vertices =  vert_df)

set.seed(11)
ggraph(ddit3_igraph, layout = 'fr') + 
  geom_edge_link(aes(edge_color = mode_of_reg, edge_alpha = guidance_ident)) +
  geom_node_point(aes(color = identity, alpha = guidance_ident)) +
  geom_node_label(aes(filter = guidance_ident == 'Guidance Mol', label = name,
                      color = identity), size = 3) +
  scale_alpha_manual(values = c('Guidance Mol' =1, 'Other' = 0.1)) +
  scale_edge_alpha_manual(values = c('Guidance Mol' =1, 'Other' = 0.1)) +
  scale_edge_color_manual(values = scales::hue_pal()(2) %>% set_names(c('Positive','Negative'))) +
  scale_color_manual(values = scales::hue_pal()(2) %>% set_names(c('High Stress','Low Stress')),
                     na.value = 'black') +
  guides(edge_alpha = 'none', alpha = 'none',
         edge_color = guide_legend(title = 'Mode of Reg'),
         color = guide_legend(title = 'Identity')) +
  theme_void()
