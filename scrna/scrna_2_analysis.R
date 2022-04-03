#Implement Single Cell Analysis

library(Seurat)
library(tidyverse)
library(reshape2)
library(Matrix)
library(patchwork)

#Load in Files
zonal_annotation <- read.table('ORs-by-zone.txt', sep = '\t', header = F,
                               col.names = c('gene_name', 'zone'))
iRFP_DE_High_vs_Low <- read_tsv('DE_iRFP_Bright_vs_Dim.tsv')

#Included here to show how files_dir was generated. Use the .rds file included, below.
#files_dir <- c('/path/to/10x/dir/1', '/path/to/10x/dir/1') %>%
#  map(.f = function(filepath_) {
#    extracted_name <- str_extract(filepath_, '(Tex15|Lhx2)?WT|wt_rep[0-9]|gg8tta')
#    Read10X(data.dir = filepath_) %>%
#      CreateSeuratObject(project = extracted_name, min.cells = 3, min.features = 200) %>%
#      list() %>% set_names(extracted_name)
#  }) %>% unlist(recursive = F) 
files_dir <- readRDS('raw_scrnaseq_seurat_ls.rds') #Raw scRNA-seq counts from Lomv Lab Datasets

axon_guidance_go_terms <- read_tsv('axon_GOs.tsv', 
                                   col_names = c('go_id','go_description')) #from Amigo2 search for "axon+guidance", see methods.

activity_dep <- map2(
  .x = c('upregulated_activity_wang.tsv',
         'downregualted_activity_wang.tsv'),
  .y = c('upregulated_activity', 'downregulated_activity'),
  .f = function(.x,.y) {
    scan(.x, sep = '\t', what = 'chr') %>% 
      matrix(nrow = 3) %>% t() %>% 
      .[,1] %>% .[. != ''] %>%
      list() %>% set_names(.y)
  }) %>% unlist(recursive = F)

#Integrate the scRNAseq Data using SCTransform Pipeline 
#see also https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
files_dir <- lapply(X = files_dir, FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = files_dir, nfeatures = 3000)
files_dir <- PrepSCTIntegration(object.list = files_dir, anchor.features = features)

anchors_ <- FindIntegrationAnchors(object.list = files_dir, normalization.method = "SCT",
                                   anchor.features = features)

all_wt_cells <- IntegrateData(anchorset = anchors_, normalization.method = "SCT")
remove(files_dir, anchors_) #save RAM.

#Run Scaling/PCA/UMAP/Clustering on Integrated Data
DefaultAssay(all_wt_cells) <- "integrated"
all_wt_cells <- all_wt_cells %>%
  # ScaleData() %>% #it is VERY clear that Seurat v3.1 you do NOT run ScaleData. see https://satijalab.org/seurat/archive/v3.1/integration.html, SCTransform tab
  RunPCA() %>%
  RunUMAP(reduction = 'pca', dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5, dims = 1:30)

#Look at Clustering and Confirm Batch Effects Regressed Out
DimPlot(all_wt_cells, label = T) + theme_bw() + theme(legend.position = 'none')
DimPlot(all_wt_cells, label = T, split.by = 'orig.ident')
DimPlot(all_wt_cells, group.by = 'orig.ident') #solved the batch effect problem.

#Set Default Assay: Will use Use SCT for Visualization and Most Other purposes
DefaultAssay(all_wt_cells) <- "SCT"

#Show Atf5 exprs rel. to INP,iOSN,mOSN markers
{map(c('Atf5','Neurog1','Gng8','Omp'), .f = function(x) {
  FeaturePlot(all_wt_cells, features = x) +
    theme_bw() +
    theme(legend.position = 'bottom') + guides(color = guide_colorbar(direction = 'horizontal'))
}) %>% purrr::reduce(`+`)} + plot_layout(nrow = 1)

#Find mOSNs
FeaturePlot(all_wt_cells, 'Omp', label = T) + theme_bw() +
  theme(plot.title = element_blank(), legend.position = 'bottom') +
  labs(color = 'Omp') + guides(color = guide_colorbar(direction = 'horizontal'))

#Delineate Omp Cells
omp_clusters <- c('0','3','4','5','7')
omp_cells <- all_wt_cells@meta.data %>% 
  rownames_to_column(var = 'cell_id') %>% 
  filter(seurat_clusters %in% omp_clusters) %>% pull(cell_id) 

#Determine chosen OR in each mOSN
#Approach examines raw counts for top OR and 2nd-highest OR in each cell. 
#Calls "chosen OR" if top counts/2nd > min_OR_enrichment
min_OR_enrichment <- 2 
or_id_calls <- all_wt_cells@assays$SCT@counts[,omp_cells] %>% 
  .[str_detect(rownames(.),'Olfr'),] %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  melt(id.vars = c('gene_name'), variable.name = 'cell_id') %>%
  filter(value != 0) %>% 
  group_by(cell_id) %>%
  summarize('chosen_OR' = gene_name[which.max(value)],
            'chosen_OR_cts' = max(value),
            'total_OR_cts' = sum(value),
            'total_ORs_detected' = length(value),
            'frac_top_OR' = chosen_OR_cts/total_OR_cts,
            'frac_second_OR' = sort(value, decreasing = T)[2]/total_OR_cts) %>%
  mutate('identity' = ifelse(frac_top_OR/frac_second_OR > min_OR_enrichment | 
                               total_ORs_detected == 1, 'Consensus Reached','No Consensus'),
         'frac_second_OR' = replace_na(frac_second_OR,0)) #NAs are from cells with only one OR

#Add chosen OR annotations (+stress level&zone of chosen OR) back to meta.data
all_wt_cells@meta.data <- all_wt_cells@meta.data %>%
  rownames_to_column(var = 'cell_id') %>%
  left_join(or_id_calls[or_id_calls$identity == 'Consensus Reached',
                        c('cell_id','chosen_OR')], 
            by = 'cell_id') %>% #only annotate cells where consensus could be reached
  left_join(iRFP_DE_High_vs_Low[,c('gene_name','log2FoldChange','significance')], by = c('chosen_OR' = 'gene_name')) %>%
  mutate('iRFP_ID' = ifelse(log2FoldChange > 0, 'iRFP_Bright','iRFP_Dim')) %>%
  left_join(zonal_annotation, by = c('chosen_OR' = 'gene_name')) %>%
  column_to_rownames(var = 'cell_id')

#Create new Seurat object w/ ONLY mOSNs that could be annotated stress level/zone
omp_only <- all_wt_cells[,!(is.na(all_wt_cells@meta.data$iRFP_ID) | is.na(all_wt_cells@meta.data$zone))]
#omp_only <- readRDS('mOSNs_only_Seuratobj.rds') #included on Zenodo as convenience.

#mOSN UMAP by Bright (stress score >0) vs. Dim (<0) Grouping
DimPlot(omp_only, group.by = 'iRFP_ID') + 
  theme_bw() +
  coord_cartesian(xlim = c(-2.5,11), ylim = c(-7.5,6.5)) +
  theme(legend.position = 'bottom') +
  guides('color' = guide_legend(nrow = 2, override.aes = list(size = 4))) +
  ggtitle('Stress Level')

#mOSN UMAP by stress score rank
irfp_ranks <- iRFP_DE_High_vs_Low %>% 
  filter(str_detect(gene_name,'Olfr')) %>%
  mutate('stress_score_rank' = rank(log2FoldChange)) %>% 
  dplyr::select(gene_name,stress_score_rank)

ggplot(data = omp_only@reductions$umap@cell.embeddings %>%
         as.data.frame() %>%
         rownames_to_column(var = 'cell_id') %>%
         left_join(rownames_to_column(omp_only@meta.data, var = 'cell_id'), 
                   by = 'cell_id') %>%
         left_join(irfp_ranks, by = c('chosen_OR' = 'gene_name')),
       aes(x = UMAP_1, y = UMAP_2, color = stress_score_rank)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient2(low = scales::hue_pal()(2)[2], mid = 'lightgrey', high = scales::hue_pal()(2)[1], 
                        midpoint = max(irfp_ranks$stress_score_rank)/2) +
  coord_cartesian(xlim = c(-2.5,11), ylim = c(-7.5,6.5)) +
  theme(legend.position = 'bottom') +
  guides('color' = guide_colorbar(direction = 'horizontal')) +
  ggtitle('Stress Level')

#mOSN UMAP by Zonal Identity
change_zonal_annotation <- function(x) {
  new <- x
  new@meta.data$zone <- ifelse(new@meta.data$zone == 'fishOR','Class I',as.character(new@meta.data$zone))
  new
} #helper function to quickly swap fishOR (annotation label) -> Class I (used in paper)

DimPlot(change_zonal_annotation(omp_only), group.by = 'zone') + 
  theme_bw() +
  coord_cartesian(xlim = c(-2.5,11), ylim = c(-7.5,6.5)) +
  theme(legend.position = 'bottom') +
  guides('color' = guide_legend(nrow = 2, override.aes = list(size = 4))) +
  ggtitle('Zone') +
  scale_color_manual(values = c(scales::hue_pal()(9),'black') %>% 
                       set_names(change_zonal_annotation(omp_only)@meta.data$zone %>%
                                   factor() %>% levels()))

#Scran Pairwise DE analysis: find differential bright/dim genes within zone.
#Note: it is IDENTICAL to use SCT@data (corrected, logtransformmed) vs RNA@data (uncorrected). Get exact SAME gene set, almost exact same p values
scran_p_vals <- scran::pairwiseWilcox(x = omp_only@assays$SCT@data, 
                                      groups = omp_only@meta.data$iRFP_ID,
                                      block = factor(omp_only@meta.data$zone),
                                      BPPARAM = BiocParallel::MulticoreParam(10)) %>%
  .$statistics %>% .[[1]] %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_name')

#Compute Within-Zone Log2FCs Genome-Wide
#1) Log2FC Bright/Dim in Each Zone
logfcs <- map_dfr(levels(omp_only@meta.data$zone), .f = function(x) {
  bright <- rowMeans(omp_only@assays$SCT@counts[,omp_only@meta.data$zone == x & omp_only@meta.data$iRFP_ID == 'iRFP_Bright'])
  dim <- rowMeans(omp_only@assays$SCT@counts[,omp_only@meta.data$zone == x & omp_only@meta.data$iRFP_ID == 'iRFP_Dim'])
  data.frame('gene_name' = rownames(omp_only),
             'log2FoldChange' = log2((bright+1)/(dim+1)),
             'zone' = rep(x, dim(omp_only)[1]))
}) #Seurat uses log10FC, I will use log2 here. 
#Checked the rest of the formula is ok by using FindMarkers Bright/Dim on the whole thing and comparing that to the result here using log10 transform

#2) Zone Sizes -> Zonal Weights
zone_sizes <- omp_only@meta.data %>% 
  group_by(iRFP_ID, zone) %>% 
  summarize(n=n()) %>% 
  group_by(zone) %>% 
  summarize('weight' = n[1]*n[2]) %>% #weight within block = Nx*Ny (NOT '+', see scran::PairwiseWilcox)
  ungroup()

#3) Final Log2FCs, Collect P values
gene_means <- data.frame('gene_name' = rownames(omp_only),
                         'baseMean' = log10(rowMeans(omp_only@assays$SCT@counts) + 1))

avg_logfcs <- logfcs %>% 
  left_join(zone_sizes, by = 'zone') %>%
  group_by(gene_name) %>%
  summarize('log2FoldChange' = weighted.mean(log2FoldChange,weight)) %>%
  ungroup() %>%
  left_join(gene_means, by = 'gene_name') %>%
  left_join(scran_p_vals, by = 'gene_name')

#Add axon guidance annotations to DE dataframe
ensembl <-  biomaRt::useEnsembl(biomart = "ensembl",
                                dataset = "mmusculus_gene_ensembl",
                                version = "97")
axon_guidance_go <- axon_guidance_go_terms %>%
  inner_join(biomaRt::getBM(attributes = c('go_id','external_gene_name'), 
                            filters = 'go_parent_term',
                            values = .$go_id,
                            mart = ensembl,
                            useCache = F),
             by = 'go_id') #just gets actual GO_IDs of interest b/c getBM returns all go-ids associated with Genes with the given GOid

avg_logfcs <- avg_logfcs %>% 
  mutate('significance' = ifelse(FDR < 0.05, 'Significant','Not Significant'),
         'guidance_ident' = ifelse(gene_name %in% axon_guidance_go$external_gene_name |
                                     str_detect(gene_name, 'Efna5|Epha5|Pcdh(10|17|19)|Kirrel(2|3)|Cntn4|Nrp|Tenm(2|4)'),
                                   'Axon Guidance Molecule','Other'))
#avg_logfcs <- read_tsv('DE_scRNAseq_Bright_Dim_withinzone.tsv') #added on Zenodo as convenience

#Curate Set of Guidance Molecules that are DE in Bright/Dim Cells
guidance_network_toplot <- avg_logfcs %>%
  filter(guidance_ident == 'Axon Guidance Molecule' & significance == 'Significant' & 
           !(gene_name %in% c('Tubb2b','Evl','Lhx2','Foxp1', 'Rpl24','Vasp','Arhgap35', 'Isl2',
                              'Etv4'))) %>% #manually remove some genes in axon guidance GO category, but TFs/other
  mutate('stress_ident' = ifelse(log2FoldChange>0,'High Stress','Low Stress'),
         'stress_ident' = fct_relevel(stress_ident, 'High Stress'),
         'identity' = map_chr(gene_name, .f = function(x) {
           if (x %in% activity_dep$upregulated_activity) {'Upregulated Activity'}
           else if (x %in% activity_dep$downregulated_activity) {'Downregulated Activity'}
           else {'Other'}
         }),
         'identity' = as.character(fct_recode(identity,
                                              'black' = 'Other', 
                                              'blue' = 'Upregulated Activity',
                                              'red' = 'Downregulated Activity'))) %>%
  
  arrange(log2FoldChange)
#NOTE: the activity annotations here are based ONLY on Wang et al., 2017
#I manually added the annotations from Tsukahara et al., 2021 in illustrator, since there were just a few additional downreg. activity molecules
#guidance_network_toplot <- read_tsv('guidance_network.tsv')

#Make TilePlot for Axon Guidance Molecules
n_low_stress <- length(which(guidance_network_toplot$stress_ident == 'Low Stress'))
groupings_ <- rep(seq(1,nrow(guidance_network_toplot)), each = n_low_stress)
groupings_ <- groupings_[1:nrow(guidance_network_toplot)]
groupings_[(length(groupings_) - n_low_stress +1):length(groupings_)] <- max(groupings_)
groupings_ <- rev(groupings_)

plot_ls <- guidance_network_toplot %>%
  split(groupings_) %>%
  map(.f = function(x) {
    limits_ <- x$gene_name
    if (length(limits_) < n_low_stress) {
      limits_ <- c(c(rep('', n_low_stress - length(limits_))), limits_)
    }
    ggplot(data = x,
           aes(x='a',y = gene_name, fill = log2FoldChange)) +
      geom_tile(size = 1, color = 'black') +
      theme_bw() +
      coord_equal() +
      scale_fill_gradient2(low = scales::hue_pal()(2)[2], mid = 'lightgrey', high = scales::hue_pal()(2)[1],
                           midpoint = 0, limits = range(guidance_network_toplot$log2FoldChange)) +
      theme(axis.text.y = element_text(color = x$identity),
            axis.text.x = element_blank(), axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank(), panel.grid = element_blank(), 
            panel.border = element_blank(), plot.margin = margin(0,11/2,0,0,'pt')) +
      scale_y_discrete(position = 'right', limits = limits_) +
      guides('fill' = guide_colorbar(direction = 'horizontal'))
  })

(purrr::reduce(plot_ls,`+`) + plot_layout(nrow = 1))/guide_area() + plot_layout(guides = 'collect', heights = c(0.9,0.1))
#again, see above note that activity part of this only is based on Wang et al., 2017
#Had to add a couple of downreg. activity molecules manually from Tsukahara et al., 2021 after making this figure. See paper.

#Volcano Plot scRNA iRFP Bright/Dim, showing guidance molecules
ggplot(data = avg_logfcs %>%
         mutate('guidance_ident' = ifelse(gene_name %in% guidance_network_toplot$gene_name,
                                          'Guidance Molecule','Other'),
                'color_code' = ifelse(significance == 'Significant',
                                      ifelse(guidance_ident == 'Guidance Molecule',
                                             'Significant Guidance Molecule',
                                             'Significant Other'),
                                      'Not Significant'),
                'color_code' = fct_relevel(color_code, 'Not Significant','Significant Other')) %>% #ensure guidance molecules plotted last.
         arrange(color_code),
       aes(x = log2FoldChange, y = -log10(FDR), color = color_code, alpha = color_code)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('Not Significant' = 'black', 
                                'Significant Guidance Molecule' = 'blue',
                                'Significant Other' = 'red')) +
  scale_alpha_manual(values = c('Not Significant' = 0.2, 
                                'Significant Guidance Molecule' = 1,
                                'Significant Other' = 0.2)) +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  guides('color'= guide_legend(nrow = 2, override.aes = list(alpha = 1))) +
  xlab('Avg log2FC High/Low Stress Within Zone')

#Working here.
#saved image as /media/storageA/hani/WorkingDir_GithubZenodoUpload/working32922.RData
#basically the code should be set from now. Just continue saving useful intermediates
#ultiamtely will just need to change the top parts of the code to import right files from Zenodo.
#Note guidance network to plot has to be c/w previous guidance_network from bulk analysis. Does NOT include the annotations from Bob
#Will just make a note that this was added for paper, not worth re-saving everything. 





#Identify the Zonal DE molecules blocking across iRFP ID
scran_zonal_raw_data <- scran::pairwiseWilcox(x = omp_only@assays$SCT@data, 
                                              groups = omp_only@meta.data$zone,
                                              block = factor(omp_only@meta.data$iRFP_ID),
                                              BPPARAM = MulticoreParam(10)) %>%
                                              {map2_dfr(.$statistics, 
                                                        base::split(as.data.frame(.$pairs), 
                                                                    seq(1, nrow(.$pairs))), 
                                                        .f = function(a,b) {
                                                          as.data.frame(a) %>%
                                                            rownames_to_column(var = 'gene_name') %>%
                                                            mutate('first' = b$first,
                                                                   'second' = b$second)
                                                        })} 

zonal_log2fcs_by_irfp <- map_dfr(unique(omp_only@meta.data$iRFP_ID), .f = function(irfp_lev) {
  print(paste0('iRFP ID:',irfp_lev))
  zonal_comps <- scran_zonal_raw_data %>% 
    filter(first == '1') %>%
    group_by(first, second) %>% summarize()
  
  map2_dfr(zonal_comps$first, zonal_comps$second, .f = function(a, b) {
    print(paste0('Comparison:',a,'.',b))
    a_means <- rowMeans(omp_only@assays$SCT@counts[,omp_only@meta.data$zone == a & omp_only@meta.data$iRFP_ID == irfp_lev])
    b_means <- rowMeans(omp_only@assays$SCT@counts[,omp_only@meta.data$zone == b & omp_only@meta.data$iRFP_ID == irfp_lev])
    data.frame('gene_name' = rownames(omp_only),
               'log2FoldChange' = log2((a_means+1)/(b_means+1)),
               'first' = rep(a,dim(omp_only)[1]),
               'second' = rep(b,dim(omp_only)[1]),
               'iRFP_ID' = rep(irfp_lev, dim(omp_only)[1]))
  })
})

cell_count_by_zone_id <- omp_only@meta.data %>%
  group_by(zone, iRFP_ID) %>%
  summarize(n=n())

zonal_scran_final_out <- zonal_log2fcs_by_irfp %>%
  left_join(cell_count_by_zone_id, by = c('first' = 'zone','iRFP_ID')) %>%
  dplyr::rename('n_first' = 'n') %>%
  left_join(cell_count_by_zone_id, by = c('second' = 'zone', 'iRFP_ID')) %>%
  dplyr::rename('n_second' = 'n') %>%
  mutate('weight' = n_first * n_second) %>% #see scran::pairwiseWilcox, weight for block should be Nx*Ny NOT '+'
  group_by(gene_name, first, second) %>%
  summarize('log2FoldChange' = weighted.mean(log2FoldChange, w = weight)) %>%
  left_join(scran_zonal_raw_data, by = c('gene_name','first','second')) %>%
  left_join(gene_means, by = 'gene_name')

#Make Plots
er_genes <- read_tsv('/media/storageA/hani/11-26_RNA_p7_Acetophenone_iRFPBrightDim/er_go_terms.tsv', 
                     col_names = c('go_id','go_description')) %>%
  filter(str_detect(go_description,'chaperone')) %>%
  inner_join(biomaRt::getBM(attributes = c('go_id','external_gene_name'), 
                            filters = 'go',
                            values = .$go_id,
                            mart = ensembl),
             by = 'go_id')

#For Supplemental Figure, Show Zone 5 and FishORs
ggplot(data = zonal_scran_final_out %>%
         filter(second %in% c('5','fishOR')) %>%
         mutate('identity' = ifelse(gene_name %in% er_genes$external_gene_name,
                                    'Chaperone','Other'),
                'sig' = ifelse(FDR < 0.05, 'Significant','Not Significant'),
                'log2FoldChange' = -log2FoldChange, #it is right now z1/z5 etc. Want Z5/Z1
                'second' = fct_recode(second, 'Class I' = 'fishOR',
                                      'Zone 5' = '5'),
                'second' = fct_relevel(second, 'Class I','Zone 5')),
       aes(x = log2FoldChange, y = -log10(FDR), color = sig, alpha = identity)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_point() +
  facet_grid(cols = vars(identity), rows = vars(second)) +
  ggrepel::geom_label_repel(data = . %>%
                              filter(sig == 'Significant' & identity == 'Chaperone') %>%
                              group_by(second) %>%
                              arrange(FDR, .by_group = T) %>%
                              dplyr::slice(1:10) %>%
                              ungroup(),
                            aes(label = gene_name), show.legend = F, min.segment.length = 0,
                            size = 3, box.padding = 1) +
  theme_bw() +
  scale_color_manual(values = c('Significant' = 'red','Not Significant' = 'black')) + 
  scale_alpha_manual(values = c('Chaperone' = 1,'Other' = 0.2)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.title = element_blank()) + 
  xlab('Mean Expression') + ylab('Avg Log2FC vs. Zone 1')

#For curiosity, same thing with axon guidance molecules
ggplot(data = zonal_scran_final_out %>%
         filter(second %in% c('5','fishOR')) %>%
         mutate('identity' = ifelse(gene_name %in% axon_guidance_go$external_gene_name |
                                      str_detect(gene_name, 'Efna5|Epha5|Pcdh(10|17|19)|Kirrel(2|3)|Cntn4|Nrp|Tenm(2|4)'),
                                    'Guidance Molecule','Other'),
                'sig' = ifelse(FDR < 0.05, 'Significant','Not Significant'),
                'log2FoldChange' = -log2FoldChange, #it is right now z1/z5 etc. Want Z5/Z1
                'second' = fct_recode(second, 'Class I' = 'fishOR',
                                      'Zone 5' = '5'),
                'second' = fct_relevel(second, 'Class I','Zone 5')),
       aes(x = log2FoldChange, y = -log10(FDR), color = sig, alpha = identity)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_point() +
  facet_grid(cols = vars(identity), rows = vars(second)) +
  ggrepel::geom_label_repel(data = . %>%
                              filter(sig == 'Significant' & identity == 'Guidance Molecule') %>%
                              group_by(second) %>%
                              arrange(FDR, .by_group = T) %>%
                              dplyr::slice(1:10) %>%
                              ungroup(),
                            aes(label = gene_name), show.legend = F, min.segment.length = 0) +
  theme_bw() +
  scale_color_manual(values = c('Significant' = 'red','Not Significant' = 'black')) + 
  scale_alpha_manual(values = c('Guidance Molecule' = 1,'Other' = 0.2)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.title = element_blank()) + 
  xlab('Mean Expression') + ylab('Avg LogFC vs. Zone 1')
