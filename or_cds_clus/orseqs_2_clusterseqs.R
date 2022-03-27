#Cluster the OR Protien Sequences 

library(DECIPHER)
library(tidyverse)
library(reshape2)
library(ggdendro)
library(patchwork)

#Load files
setwd('./')
parsed_fasta <- readDNAStringSet(filepath = 'OR_CDS.fa')
iRFP_levels <- read_tsv('DE_iRFPBright_vs_Dim.tsv')
zones <- read_tsv('ORs-by-zone.txt', col_names = c('gene_name','zone'))

#Translate Sequences in proper frame
protein_seqs<- translate(parsed_fasta[(width(parsed_fasta)/3)%%1==0]) #seq must be divisible by 3...

#Run the protein alignment with DECIPHER
aa_aligned<- AlignSeqs(protein_seqs)
aa_staggered <- StaggerAlignment(aa_aligned)
aa_dist <- DistanceMatrix(aa_staggered, type = 'dist')

#Remove Na's from Distance Matrix and Convert
dist_no_nas <- as.matrix(aa_dist)
dist_no_nas <- dist_no_nas[rowSums(is.na(dist_no_nas)) == 0, colSums(is.na(dist_no_nas)) == 0, drop = FALSE]
dist_no_nas <- as.dist(dist_no_nas)

#Run MDS
mds_out <- cmdscale(dist_no_nas, k = 30) 
set.seed(2)
kmeans_out <- kmeans(mds_out, centers = 6, nstart = 20)

mds_to_plot <- mds_out %>% 
  as.data.frame() %>%
  set_names(paste0('PC_', seq(1,ncol(mds_out)))) %>%
  rownames_to_column(var = 'or_name') %>%
  mutate('gene_name' = str_extract(or_name,'Olfr[0-9]+(\\-ps([0-9]+)?)?'),
         'cluster' = kmeans_out$cluster) %>%
  left_join(dplyr::select(iRFP_levels,gene_name,log2FoldChange),
            by = 'gene_name') %>%
  mutate('iRFP_ID' = ifelse(log2FoldChange > 0, 'iRFP_Bright','iRFP_Dim')) %>%
  left_join(zones, by = 'gene_name')
#write_tsv(mds_to_plot, 'mds_irfp_zone_orseqs.tsv')

#Plot, mapping stress score rank to color
ggplot(data = mds_to_plot, 
       aes(x = PC_1, y = PC_2, color = rank(log2FoldChange))) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_color_gradient2(low = scales::hue_pal()(2)[2], mid = 'lightgrey', high = scales::hue_pal()(2)[1], 
                        midpoint = nrow(mds_to_plot)/2) +
  theme(legend.position = 'bottom', legend.title = element_text(vjust = 1)) +
  guides('color' = guide_colorbar(title = 'Stress Score Rank', label = F, ticks = F,
                                  barwidth = 3, barheight = 0.5))

#Plot, mapping zone to color
ggplot(data = mds_to_plot %>%
         filter(!is.na(zone)) %>%
         mutate('class' = ifelse(zone == 'fishOR','Class I', 'Class II'),
                'Zone' = ifelse(class == 'Class I','Class I',zone)), 
       aes(x = PC_1, y = PC_2, color = Zone)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  guides('color' = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  scale_color_manual(values = c(scales::hue_pal()(9),'black') %>% 
                       set_names(c(sort(unique((mds_to_plot$zone)))[1:9],'Class I')))
