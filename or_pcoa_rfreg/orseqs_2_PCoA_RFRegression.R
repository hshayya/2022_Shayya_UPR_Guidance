#PCoA Visualization and RF Regression of OR AA Sequences

#############################################################
####                  Setup & Imports                    ####
#############################################################

#Set up the Environment 

library(DECIPHER)
library(tidyverse)
library(reshape2)
library(ggdendro)
library(patchwork)
library(randomForest)
library(caret)
library(doParallel)

#Load FASTA and Annotations for Use Later
parsed_fasta <- readDNAStringSet(filepath = 'OR_CDS.fa')

iRFP_levels <- read_tsv('DE_iRFPBright_vs_Dim.tsv')
stress_scores <- iRFP_levels %>% 
  filter(str_detect(gene_name,'Olfr')) %>%
  mutate('identity' = ifelse(log2FoldChange > 0,'High Stress','Low Stress'),
         'stress_score_rank' = rank(log2FoldChange))

zones <- read_tsv('ORs-by-zone.txt', 
                  col_names = c('gene_name','zone'))

source('snakeplotter_mod.R') #have to make sure this uploaded Github

#M71/M72 Mutant and B6/129 Phenos from Feinstein and Mombaerts 2004
momb_phenos <- c('momb_gfp_lacz_phenos.csv',
                 'momb_lacz_phenos.csv',
                 'momb_strain_phenos.csv') %>%
  map_dfr(.f = function(x) {
    read_csv(x) %>% melt(id.vars = 'X1', variable.name = 'prey') %>%
      dplyr::rename(bait = X1)
  }) %>%
  dplyr::filter(!is.na(value)) %>%
  mutate('comp' = map2_chr(bait, prey, .f = function(a, b) {paste0(sort(c(as.character(a),
                                                                          as.character(b))), 
                                                                   collapse = '')})) %>%
  group_by(comp) %>%
  summarize('value' = unique(value))

#############################################################
####               PCoA & Visualization                  ####
#############################################################

#Translate Sequences in proper frame
protein_seqs<- translate(parsed_fasta[(width(parsed_fasta)/3)%%1==0])

#Run the protein alignment with DECIPHER
aa_aligned<- AlignSeqs(protein_seqs)
aa_staggered <- StaggerAlignment(aa_aligned)
aa_dist <- DistanceMatrix(aa_staggered, type = 'dist')
aa_dendro <- IdClusters(aa_dist, type = 'dendrogram')

#Cluster the Alignment 
dist_no_nas <- as.matrix(aa_dist)
dist_no_nas <- dist_no_nas[rowSums(is.na(dist_no_nas)) == 0, colSums(is.na(dist_no_nas)) == 0, drop = FALSE]
dist_no_nas <- as.dist(dist_no_nas)

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

#Re-name clusters by stress score, makes it easier for visualization
cluster_remapping <- mds_to_plot %>% 
  group_by(cluster) %>%
  summarize('avg_score' = mean(log2FoldChange, na.rm = T)) %>%
  arrange(avg_score) %>%
  mutate('new_cluster' = seq(1,nrow(.)))

mds_to_plot <- mds_to_plot %>%
  left_join(dplyr::select(cluster_remapping, cluster, new_cluster))

#write_tsv(mds_to_plot, 'OR_PCs_stressscores_zone_cluster.tsv')

#PCoA Plot: color ~ Cluster IDs
ggplot(data = mds_to_plot, 
       aes(x = PC_1, y = PC_2, color = as.character(new_cluster))) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  labs(color = 'Cluster') +
  theme(legend.position = 'bottom') + guides(color = guide_legend(nrow = 1))

#PCoA Plot: color ~ stress score rank
ggplot(data = mds_to_plot, 
       aes(x = PC_1, y = PC_2, color = rank(log2FoldChange))) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_color_gradient2(low = scales::hue_pal()(2)[2], mid = 'lightgrey', high = scales::hue_pal()(2)[1], 
                        midpoint = nrow(mds_to_plot)/2) +
  theme(legend.position = 'bottom', legend.title = element_text(vjust = 1)) +
  guides('color' = guide_colorbar(title = 'Stress Score Rank', label = F, ticks = F,
                                  barwidth = 3, barheight = 0.5))

#Vln Plot stress score distribution ~ cluster
ggplot(data = mds_to_plot, 
       aes(x = as.character(new_cluster), y = log2FoldChange, 
           fill = as.character(new_cluster))) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_violin() +
  geom_boxplot(width = 0.2, fill = 'white') +
  theme_bw() +
  theme(legend.position = 'none') +
  xlab('Cluster') + ylab('Stress Score')

#PCoA Plot: color ~ zone
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

#Collapse Annotations
annotation_df <- enframe(kmeans_out$cluster, name = 'seq', value = 'cluster') %>%
  mutate('gene_name' = str_extract(seq, '.*?(?=_)')) %>%
  left_join(dplyr::select(cluster_remapping, cluster, new_cluster)) %>%
  left_join(zones, by = 'gene_name') %>%
  left_join(dplyr::select(stress_scores, gene_name, identity, stress_score = log2FoldChange))

#write_tsv(annotation_df, 'ORseqs_annotations.tsv')

#############################################################
####             Random Forest Regression                ####
#############################################################

#We will need to add Mombaerts Sequences to MSA if we want to ultimately predict stress scores for them w/ RF Model
#All of the sequence information below is from Feinstein and Mombaerts 2004 Cell Paper, Methods & Figures

#Start w/ M72 B6 sequence
m72_seq <- as.vector(as.matrix(protein_seqs['Olfr160_ENSMUSG00000061165_ENSMUST00000215727']))

#Function to adjust sequences
adjust_protein_seq <- function(input_vec, pos, aa_string) {
  new_vec <- input_vec
  aa_vec <- str_split(aa_string,'')[[1]]
  new_vec[pos] <- aa_vec
  new_vec
}

#M72 129 Sequence
m72_129_seq <- adjust_protein_seq(m72_seq, c(37,284), 'VF')

#M71 129 sequence (M71 B6 = pseudogene). 11 mutations from M72 129 Sequence.
m71_seq <- adjust_protein_seq(m72_129_seq, c(2,15,28,36,145,151,157,161,194,205,307), 'TGFVAVAGVDK')

#ORs w/ 129/B6 Polymorphisms. Below are the B6 sequences
p2_seq <- as.vector(as.matrix(protein_seqs['Olfr17_ENSMUSG00000073897_CUFFORT_1594']))
m50_seq <- as.vector(as.matrix(protein_seqs['Olfr6_ENSMUSG00000036647_CUFFORT_1781']))
p4_seq <- as.vector(as.matrix(protein_seqs['Olfr714_ENSMUSG00000049674_ENSMUST00000214429']))

#New Sequences from Feinstein and Mombaerts 2004
new_seqs <-
  list('Olfr151_129' = m71_seq,
       'MutA_Momb' = adjust_protein_seq(m71_seq, c(2,15,28,36), 'ARLI'),
       'MutB_Momb' = adjust_protein_seq(m71_seq, c(145,151,157,161,194,205), 'VMTSIN'),
       'MutC_Momb' = adjust_protein_seq(m71_seq, c(307), 'R'),
       'MutD_Momb' = adjust_protein_seq(m71_seq, c(145,151,157,161), 'VMTS'),
       'MutE_Momb' = adjust_protein_seq(m71_seq, c(145,151), 'VM'),
       'MutF_Momb' = adjust_protein_seq(m71_seq, c(2,15,145,151), 'ARVM'),
       'MutG_Momb' = adjust_protein_seq(m71_seq, c(15,151), 'RM'),
       'MutH_Momb' = adjust_protein_seq(m71_seq, c(157,161), 'TS'),
       'MutI_Momb' = adjust_protein_seq(m71_seq, c(157), 'T'),
       'MutJ_Momb' = adjust_protein_seq(m71_seq, c(205), 'N'),
       'P2_129' = adjust_protein_seq(p2_seq, 282, 'I'),
       'M50_129' = adjust_protein_seq(m50_seq, c(86,161), 'RI'),
       'P4_129' = adjust_protein_seq(p4_seq, c(142), 'H'),
       'M72_129' = m72_129_seq) %>% 
  imap_dfr(.f = function(x,y) {
    data.frame('position_unaligned' = seq(1,length(x)),
               'value' = x,
               'new_seq' = y)
  })

#Run the re-alginment
protein_seqs_wmombaerts <- new_seqs %>% 
  group_by(new_seq) %>% 
  group_map(function(.x,.y) {paste0(.x$value, collapse='') %>% 
      set_names(unlist(.y))}) %>% 
  unlist(recursive = F) %>%
  AAStringSet() %>%
  {append(protein_seqs,.)}

allseqs_aligned<- AlignSeqs(protein_seqs_wmombaerts)
allseqs_staggered <- StaggerAlignment(allseqs_aligned)

#Biostrings::writeXStringSet(allseqs_staggered, 'or_momb_msa.fasta') 

#Convert MSA to dataframe
allseqs_df <- as.data.frame(as.matrix(allseqs_staggered)) %>%
  rownames_to_column(var = 'sequence') %>%
  melt(id.vars = 'sequence', variable.name = 'position') %>%
  mutate('position' = as.numeric(str_extract(position, '[0-9]+')))

#Show fraction of '-' in MSA - really only want positions where most ORs have an AA...
ggplot(data = allseqs_df %>%
         group_by(position, value) %>%
         summarize(n=n()) %>%
         summarize('frac_empty' = 1-sum(n[value!='-'])/sum(n)), 
       aes(x = frac_empty)) +
  geom_histogram() +
  theme_bw()

#Non-'-' positions in the MSA
allseqs_pos_to_include <- allseqs_df %>%
  group_by(position, value) %>%
  summarize(n=n()) %>%
  summarize('frac_empty' = 1-sum(n[value!='-'])/sum(n)) %>%
  filter(frac_empty < 0.9) %>%
  ungroup() %>%
  pull(position)  

#Prep dataframe for RF by adding stress scores
allseqs_df_forcaret <- allseqs_df %>%
  filter(position %in% allseqs_pos_to_include) %>% #only include positions where at least 10% ORs have an AA
  dcast(sequence ~ position, value.var = 'value') %>%
  left_join(dplyr::select(annotation_df, seq, stress_score), by = c('sequence' = 'seq')) %>%
  rename_at(.vars = vars(-sequence, -stress_score), .funs = function(x) {paste0('V',x)}) %>%
  mutate_at(.vars = vars(-stress_score), .funs = factor) %>%
  column_to_rownames(var = 'sequence')

#Train RF Model
train_df <- allseqs_df_forcaret[!is.na(allseqs_df_forcaret$stress_score),]

params <- trainControl(method = 'repeatedcv',
                       number = 10,
                       repeats = 3, 
                       search = 'grid')

tunegrid <- expand.grid(.mtry=c(3,20,50,100))

cores <- makeCluster(10)
registerDoParallel(cores = cores)

modellist <- list()
for (ntree in c(100,500,1000)){
  set.seed(123)
  print(paste0('Working on: ', ntree, ' trees'))
  fit <- train(x = train_df[,1:322],
               y = train_df[,323], 
               method = 'rf',
               metric = 'RMSE',
               trControl = params,
               tuneGrid = tunegrid,
               ntree = ntree, 
               importance = TRUE)
  key <- toString(ntree)
  modellist[[key]] <- fit
}

stopCluster(cores)

#View Grid Search Results
ggplot(data = map_dfr(names(modellist), .f = function(x) {
  modellist[[x]]$results %>% mutate('ntree' = x)
}) %>%
  mutate('upper_bound' = RMSE + RMSESD,
         'lower_bound' = RMSE - RMSESD,
         'ntree' = fct_relevel(ntree, '100','500')),
aes(x = mtry, y = RMSE, color = ntree, group = ntree)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = 'bottom') + 
  guides(color = guide_legend(nrow = 1))

#Decide to use ntree = 500, mtry = 50. Diminishing returns after this point
model_allseqalign <- modellist$`500`
model_allseqalign <- update(model_allseqalign, param = list(mtry = 50))

#saveRDS(model_allseqalign,'rf_model.rds')

#Plot Model Fit Predicted Scores ~ Measured Scores
ggplot(data = data.frame('obs' = train_df[,323],
                         'pred' = predict(model_allseqalign)),
       aes(x = obs, y = pred)) +
  geom_point(alpha = 0.5) +
  stat_function(fun = function(x) {x}, linetype = 'dashed') +
  geom_smooth(method = 'lm', se = T, color = 'blue', fill = 'blue', alpha = 0.2) +
  theme_bw() + 
  geom_text(data = data.frame('lab' = paste0('R2=',round(model_allseqalign$results[model_allseqalign$results == 50,
                                                                                   'Rsquared'],2))), 
            aes(label = lab), 
            x = -Inf, y = Inf, hjust = 0, vjust = 1, inherit.aes = F) +
  xlab('Measured Stress Score') + ylab('RF Predicted Stress Score')

########################################################
####             Variable Importance                ####
########################################################

#Variable Importance from RF Model
var_importance <- importance(model_allseqalign$finalModel) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'position')

#Is there are relationship between Var Importance and Shannon Entropy?
shannon.entropy <- function(p)
{
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

entropy_by_pos <- train_df %>% 
  rownames_to_column(var = 'sequence') %>%
  melt(id.vars = c('sequence', 'stress_score')) %>%
  mutate('position' = as.numeric(str_extract(variable,'[0-9]+'))) %>%
  group_by(position, value) %>%
  summarize(n=n()) %>%
  mutate('frac_composition' = n/sum(n)) %>%
  summarize('entropy' = shannon.entropy(frac_composition))

ggplot(data = var_importance %>%
         mutate('position' = as.numeric(substr(position,2,nchar(position)))) %>%
         left_join(entropy_by_pos, by = 'position') %>%
         filter(position %in% allseqs_pos_to_include),
       aes(x = entropy, y = `%IncMSE`)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(method = 'lm', color = 'blue', se = T, alpha = 0.2) +
  theme(legend.title = element_blank(), legend.position = 'bottom') +
  xlab('Shannon Entropy in MSA')

#Regress Var Importance by Shannon Entropy to get Corrected %IncMSE
regressed_varimportance <- var_importance %>%
  mutate('position' = as.numeric(substr(position,2,nchar(position)))) %>%
  left_join(entropy_by_pos, by = 'position') %>%
  mutate('corrected_incMSE' = residuals(lm(formula = `%IncMSE` ~ entropy, 
                                           data = .)))

#write_tsv(regressed_varimportance,'rf_regressed_var_imp.tsv')

top_residues_regression <- regressed_varimportance %>%
  filter(position %in% allseqs_pos_to_include) %>%
  dplyr::arrange(desc(corrected_incMSE)) %>%
  dplyr::slice(1:15) %>%
  mutate('rank' = seq(1,nrow(.))) %>%
  dplyr::select(position,rank)

#Point Plot %IncMSE ~ Shannon Entropy, Top Residues by corrected %IncMSE in red
ggplot(data = var_importance %>%
         mutate('position' = as.numeric(substr(position,2,nchar(position)))) %>%
         left_join(entropy_by_pos, by = 'position') %>%
         filter(position %in% allseqs_pos_to_include) %>%
         mutate('ident' = ifelse(position %in% top_residues_regression$position[1:10],
                                 'Top Residue','Other')),
       aes(x = entropy, y = `%IncMSE`, color = ident)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(method = 'lm', color = 'blue', se = T, alpha = 0.2) +
  scale_color_manual(values = c('Top Residue' = 'red','Other' = 'black')) +
  theme(legend.title = element_blank(), legend.position = 'bottom') +
  xlab('Shannon Entropy in MSA')

##############################################
####             Snakeplots               ####
##############################################

#P2 Annotations- from Uniprot Q9JKA6
olfr17_annot <- rbind(c('Tm1',26,50),
                      c('Tm2',62,83),
                      c('Tm3',103,121),
                      c('Tm4',142,165),
                      c('Tm5',201,223),
                      c('Tm6',244,262),
                      c('Tm7',274,293)) %>%
  data.frame() %>%
  set_names(c('annot','x','xend')) %>%
  mutate_at(.vars = vars(-annot), .funs = function(x) {as.numeric(as.character(x))})

#Lookup Table MSA Positions <-> P2 Positions
lookup_table_p2 <- allseqs_df %>%
  filter(str_detect(sequence,'Olfr17_')) %>%
  arrange(position) %>%
  mutate('unaligned_position' = cumsum(value != '-'),
         'unaligned_position' = ifelse(value == '-', NA, unaligned_position)) %>%
  filter(!is.na(unaligned_position)) %>%
  dplyr::select(position, unaligned_position)

#Set up for TM/EC/IC Domain Call for P2 SnakePlot
tm_doms <- olfr17_annot$xend - olfr17_annot$x + 1
ec_doms <- (c(olfr17_annot$x, width(protein_seqs['Olfr17_ENSMUSG00000073897_CUFFORT_1594']))-1) -
  (c(0, olfr17_annot$xend)) 

for_plot <- vector(class(tm_doms), length(c(tm_doms, ec_doms)))
for_plot[c(TRUE, FALSE)] <- ec_doms
for_plot[c(FALSE, TRUE)] <- tm_doms

#Lookup Table MSA Positions & Importance Measurements <-> P2 Positions
reg_imp_p2_lookup <- lookup_table_p2 %>%
  left_join(regressed_varimportance,
            by = 'position') %>%
  dplyr::arrange(unaligned_position)

#M71 Annotations- from Uniprot Q60893
m71_annot <- rbind(c('Tm1',29,49),
                   c('Tm2',57,77),
                   c('Tm3',91,111),
                   c('Tm4',134,154),
                   c('Tm5',196,216),
                   c('Tm6',239,259),
                   c('Tm7',271,291)) %>%
  data.frame() %>%
  set_names(c('annot','x','xend')) %>%
  mutate_at(.vars = vars(-annot), .funs = function(x) {as.numeric(as.character(x))})

#Lookup Table MSA Positions <-> P2 Positions
m71_to_full_msa_lookup <- allseqs_df %>%
  filter(sequence == 'Olfr151_129' & value != '-') %>%
  dplyr::arrange(position) %>%
  mutate('unaligned_position' = seq(1, length(value)))

#Set up for TM/EC/IC Domain Call for M71 SnakePlot
tm_doms_m71 <- m71_annot$xend - m71_annot$x + 1
ec_doms_m71 <- (c(m71_annot$x, width(protein_seqs_wmombaerts['Olfr151_129']))-1) -
  (c(0, m71_annot$xend)) 

for_plot_m71 <- vector(class(tm_doms_m71), length(c(tm_doms_m71, ec_doms_m71)))
for_plot_m71[c(TRUE, FALSE)] <- ec_doms_m71
for_plot_m71[c(FALSE, TRUE)] <- tm_doms_m71

#Lookup Table MSA Positions & Importance Measurements <-> M71 Positions
reg_imp_m71_lookup <- m71_to_full_msa_lookup %>%
  left_join(regressed_varimportance,
            by = 'position') %>%
  dplyr::arrange(unaligned_position)

#Prepare Regressed P2 and M71 snakeplots, for figure
#Put red stars on the top residues, will made the labels and lines in illustrator
{do.call(snakePlot_plasma, 
         c(as.list(for_plot), #ec/tm/ic domains
           list('grey', #pcirc
                'grey', #pfill
                5, #psize
                'grey', #lcol
                paste0(ifelse(reg_imp_p2_lookup$position %in% top_residues_regression$position[1:10],
                              '*',' '),
                       collapse = '')[1:315],#as.character(protein_seqs['Olfr17_ENSMUSG00000073897_CUFFORT_1594']), #aa
                'red', #aacol
                reg_imp_p2_lookup$corrected_incMSE[1:315] #fill_indiv
           ))) +  
    guides(color = guide_colorbar(title = 'Corrected\n%IncMSE'))
  } +
  {
    do.call(snakePlot_plasma, 
            c(as.list(for_plot_m71), #ec/tm/ic domains
              list('grey', #pcirc
                   'grey', #pfill
                   5, #psize
                   'grey', #lcol
                   paste0(ifelse(reg_imp_m71_lookup$position %in% top_residues_regression$position[1:10],
                                 '*',' '),
                          collapse = '')[1:309],
                   'red', #aacol
                   reg_imp_m71_lookup$corrected_incMSE[1:309] #fill_indiv
              ))) +  
      guides(fill = guide_colorbar(title = 'Corrected\n%IncMSE'), color = 'none') +  
      theme(legend.position = 'none')
  } + plot_layout(guides = 'collect')

#Show the top residues w/ labels for P2 and M71 (added these in illustrator)

reg_imp_p2_lookup %>% 
  mutate('AA' = strsplit(as.character(protein_seqs['Olfr17_ENSMUSG00000073897_CUFFORT_1594']),
                         '')[[1]]) %>%
  right_join(top_residues_regression, by = 'position') %>%
  filter(rank <= 10) %>%
  mutate('lab' = paste0(AA,unaligned_position,'(',position,')')) %>% 
  View()

reg_imp_m71_lookup %>%
  right_join(top_residues_regression, by = 'position') %>%
  filter(rank <= 10) %>%
  mutate('lab' = paste0(value,unaligned_position,'(',position,')')) %>%
  View()

###################################################################
####          AA Composition at Top RF Positions               ####
###################################################################

#Bin ORs into deciles by stress score
annotations_decile <- annotation_df %>%
  mutate('stress_decile' = ntile(stress_score, 10))

#Plot Composition at Top 10 Positions
ggplot(data = allseqs_df %>% 
         filter(sequence %in% rownames(train_df) & 
                  position %in% top_residues_regression$position[1:10]) %>% 
         left_join(annotations_decile, by = c('sequence' = 'seq')) %>%
         group_by(position, stress_decile, value) %>%
         summarize(n=n()) %>%
         mutate('frac_composition' = n/sum(n)) %>%
         ungroup() %>%
         left_join(lookup_table_p2, by = 'position') %>%
         left_join(filter(dplyr::rename(allseqs_df, aa = value), 
                          sequence == 'Olfr17_ENSMUSG00000073897_CUFFORT_1594'), 
                   by = 'position') %>%
         mutate('label' = paste0(aa,unaligned_position,'(',position,')'),
                'position' = factor(as.character(position), levels = as.character(top_residues_regression$position[1:10]))) %>% #order by position
         arrange(position) %>% #facet var is label, so arrange by position (ordered above) to get in correct order
         mutate('label' = factor(label, levels = unique(label))), #then re-factor
       aes(x = stress_decile, y = frac_composition, fill = value)) +
  geom_col() +
  facet_wrap(facets = vars(label), nrow = 2) +
  theme_bw() + xlab('Stress Score Decile') + ylab('Fractional Composition') +
  labs(fill = guide_legend(title = 'Amino Acid')) +
  scale_x_continuous(breaks = seq(1,10,2))

#########################################################################################
####          Using RF Model to Predict M71/M72 Mutant Guidance Phenos               ####
#########################################################################################

#Predict stress scores for all mombaerts seqs
seqs_to_include <- c(as.character(unique(new_seqs$new_seq)),
                     'Olfr160_ENSMUSG00000061165_ENSMUST00000215727',
                     'Olfr17_ENSMUSG00000073897_CUFFORT_1594',
                     'Olfr6_ENSMUSG00000036647_CUFFORT_1781',
                     'Olfr714_ENSMUSG00000049674_ENSMUST00000214429')

test_df <- allseqs_df_forcaret[seqs_to_include,-323]

stress_preds <- predict(model_allseqalign, test_df)

#Order sequences by clustering distance matrix of stress predictions
levs <- hclust(dist(stress_preds)) %>% {.$labels[.$order]}

#Summarize predicted differences in stress scors & phenotypes for all Feinstein/Mombaerts Comparisons
momb_stressdiffs_phenos <- as.matrix(dist(stress_preds)) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'bait') %>%
  melt(id.vars = 'bait', variable.name = 'prey') %>%
  mutate_at(.vars = vars(bait, prey), factor, levels = levs) %>% #for ordering axes
  mutate(pair = map2_chr(bait, prey, function(a,b) {
    paste(sort(c(as.character(a),as.character(b))), collapse = '')})) %>%
  group_by(pair) %>%
  dplyr::arrange(bait, .by_group = T) %>% #top triangle only
  dplyr::slice(1) %>%
  ungroup() %>%
  left_join(dplyr::select(momb_phenos, pair = comp, config = value), 
            by = 'pair') %>%
  filter(!is.na(config)) %>%
  mutate('Configuration' = fct_recode(config, 'intermixed' = 'int',
                                      'compartmentalized' = 'comp',
                                      'adjacent/accessory' = 'adj'),
         'Configuration' = fct_relevel(Configuration,'intermixed','compartmentalized'))

#For the M71/M72 Mutant Experiments, does configuration correlate with predicted difference stress scores?
#Omit cases where bait = prey, 129/B6 experiments (concern is that could have other point mutations etc.)
ggplot(momb_stressdiffs_phenos %>% 
         filter(bait != prey & !str_detect(pair, 'Olfr160|M50|P2|P4')) %>%
         mutate('Configuration' = fct_recode(Configuration,
                                             'Int' = 'intermixed', 
                                             'Comp' = 'compartmentalized',
                                             'Adj' = 'adjacent/accessory')),
       aes(x = Configuration, y = value, color = Configuration)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(width = 0.4, fill = NA) +
  theme_bw() +
  xlab('Configuration') + ylab('Stress Difference') +
  scale_color_manual(values = c('Int' = 'black',
                                'Comp' = 'blue',
                                'Adj' = 'red')) +
  theme(legend.position = 'none')

momb_stressdiffs_phenos %>% 
  filter(bait != prey & !str_detect(pair, 'Olfr160|M50|P2|P4')) %>%
  lm(formula = value ~ Configuration) %>%
  aov() %>%
  TukeyHSD()

#########################################################################################
####           Using RF Model to Predict B6/129 Mutant Guidance Phenos               ####
#########################################################################################

#Lookup tables Common Names <-> Olfr#
lookup_commons <- data.frame('common' = c('P2','M50','P4','M72','M71'),
                             'olfr' = c('Olfr17','Olfr6','Olfr714','Olfr160','Olfr151'))


#Df of stress scores for B6 and 129 Olfrs
stress_scores_df_b6_129 <- enframe(stress_preds, name = 'sequence', value = 'stress_score') %>%
  filter(!str_detect(sequence,'Momb|Olfr151')) %>%
  mutate('common' = str_extract(sequence, '.*?(?=_)')) %>%
  left_join(lookup_commons, by = 'common') %>%
  mutate('gene_name' = ifelse(is.na(olfr), common, as.character(olfr)))

#Sequences of all ORs in the B6/129 Experiment & Mapping MSA <-> Unaligned Positions
seqs_b6_129 <- allseqs_df %>%
  filter(sequence %in% stress_scores_df_b6_129$sequence & 
           value != '-') %>%
  group_by(sequence) %>%
  dplyr::arrange(position, .by_group = T) %>%
  mutate('unaligned_position' = seq(1, length(value))) %>%
  left_join(dplyr::select(stress_scores_df_b6_129, sequence, gene_name), by = 'sequence')

#Annotations of TM Domains
tm_domains <- list('Olfr6' = c(23,43,65,85,98,118,133,153, 200,220,238,258,270,290), #Uniprot P34986 
                   'Olfr160' = c(29,49,57,77,91,111, 134,154,196,216,239,259,271,291), #used M71 uniprot Q60893 b/c M72 only has 6TMs annoted -> see methods
                   'Olfr714' = c(27,51,63,83,103,124,145,165,198,227,239,262,274,293), #Uniprot Q7TRN0
                   'Olfr17' = c(26,50,62,83,103,121,142,165,201,223,244,262,274,293) #Uniprot Q9JKA6 
)

#Create the TM domain diagrams for B6/129 sequences
b6_129_diagrams <- 
  ggplot(data = imap_dfr(tm_domains, 
                         .f = function(x, nm_) {data.frame('x' = x[c(T,F)], 
                                                           'xend' = x[c(F,T)], 
                                                           'or' = nm_)}) %>%
           mutate('or' = factor(or, levels = names(tm_domains))),
         aes(x= x, xend = xend, y = or, yend = or)) +
  geom_segment(size = 5) +
  geom_segment(data = seqs_b6_129 %>% group_by(gene_name) %>% 
                 summarize('x' = 1, 
                           'xend' = max(unaligned_position)) %>%
                 dplyr::rename(or = gene_name),
               size = 1) +
  geom_point(data = seqs_b6_129 %>% 
               group_by(gene_name, unaligned_position) %>% 
               filter(length(unique(value))>1) %>% 
               summarize(),
             aes(x = unaligned_position, y = gene_name),
             color = 'red', inherit.aes = F) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank()) +
  scale_y_discrete(labels = function(x) {lookup_commons$common[match(x, lookup_commons$olfr)]})

b6_129_diagrams

#Plot Stress Differences & phenotypes from B6/129 data
b6_129_res <- 
  ggplot(data = stress_scores_df_b6_129 %>%
           group_by(gene_name) %>% 
           summarize('Stress Difference' = abs(stress_score[1]-stress_score[2])) %>%
           left_join(momb_phenos %>% 
                       filter(str_detect(comp,'129')) %>%
                       mutate('gene_name' = str_extract(comp,'Olfr[0-9]+')),
                     by  = 'gene_name') %>%
           mutate(gene_name = factor(gene_name, levels = names(tm_domains)),
                  Configuration = fct_recode(value, 'adjacent/accessory' = 'adj',
                                             'compartmentalized' = 'comp')),
         aes(x = 1, y = gene_name, fill = `Stress Difference`)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_viridis_c() + 
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = 2, fill = Configuration), color = 'black', size = 1) +
  scale_fill_manual(values = c('adjacent/accessory' = 'red', 
                               'compartmentalized' = 'blue')) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank())

b6_129_res

#B6/129 Figure for Paper
b6_129_diagrams + b6_129_res + plot_layout(guides = 'collect', widths = c(0.8,0.2))
