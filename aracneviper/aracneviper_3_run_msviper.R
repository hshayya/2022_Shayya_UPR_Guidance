library(Seurat)
library(tidyverse)
library(reshape2)
library(Matrix) #haven't exclusively tested this
library(zinbwave)
library(patchwork)

#Setup
setwd('./')

#' Run Viper Using "Default" Parameters Recommended in Vignette
#' 
#' @description Implements the pattern of rowTtest -> standardize p values -> run ttestNull with 1000 iterations -> run viper as recommended in vignette
#' @param vst_array (matrix) containing (typically) vst-transformed count data to analyze with the msviper function
#' @param bool_x (bool vector) used to subset columns of vst_array, returning sub-matrix representing condition x (numerator)
#' @param bool_y (bool vector) as bool_x, for condition y (denominator)
#' @param regulon_obj (viper regulon) to be used for analysis
#' @param use_cores (int) number of cores to be used for analysis
#' @return (list) with two items: "viper_obj" = object returned by ms viper and "df" = dataframe converted from msviper output
run_viper_v2 <- function(vst_assay, bool_x, bool_y, regulon_obj, use_cores){
  #Subset Arrays
  print('Subsetting Arrays...')
  array_x <- vst_assay[, bool_x]
  array_y <- vst_assay[, bool_y]
  
  #Compute signatures
  print('Computing Signatures...')
  signatures_ <- viper::rowTtest(x = array_x, y = array_y)
  signatures_ <- (qnorm(signatures_$p.value/2, lower.tail = FALSE) * sign(signatures_$statistic))[, 1]
  
  #Compute null model
  print('Computing Null Model...')
  null_model <- viper::ttestNull(x= array_x, y = array_y, per = 1000,
                                 repos = TRUE, verbose = T, cores = use_cores)
  
  #Run Viper
  print('Running Viper...')
  final_out <- viper::msviper(signatures_, regulon_obj, null_model, verbose = T, cores = use_cores)
  df_obj <- summary(final_out, length(final_out$es$nes)) %>% as.data.frame() %>%
    arrange(FDR)
  
  list('viper_obj' = final_out, 'df' = df_obj)
}

#Load Regulon Object
#NB: we found that the network file can be properly loaded with the '3col' format. However, the source code fails to trim the header of this file on input
#   If the header is not trimmed, the code throws NA's and fails spectacularily
#   Trimmed version of network file is loaded here & included on Zenodo
osn_regulon <- viper::aracne2regulon(afile = 'aracne_network_out_noheader.txt',
                                     eset = '272samples_OE_expression_matrix_forArachne.tsv',
                                     format = '3col')

#Load Omp Dataset
omp_only <- readRDS('mOSNs_only_Seuratobj.rds')

#Re-SCT, get all pearson residuals. Will take a WHILE.
sc_exprs_seurat_scaled <- Seurat::SCTransform(omp_only, variable.features.n = dim(omp_only)[1])  

#Run msViper on Bright vs. Dim cells in each zone. Will take a WHILE.
viper_out_late_scaled_byzone <- 
  map(levels(factor(sc_exprs_seurat_scaled@meta.data$zone)),
      .f = function(zone_) {
        run_viper_v2(sc_exprs_seurat_scaled@assays$SCT@scale.data, 
                     bool_x = sc_exprs_seurat_scaled@meta.data$iRFP_ID == 'iRFP_Bright' & 
                       sc_exprs_seurat_scaled@meta.data$zone == zone_,
                     bool_y = sc_exprs_seurat_scaled@meta.data$iRFP_ID == 'iRFP_Dim' & 
                       sc_exprs_seurat_scaled@meta.data$zone == zone_,
                     regulon_obj = osn_regulon,
                     use_cores = 8)
      }
  )

by_zone_viper_df <- 
  map2(viper_out_late_scaled_byzone, levels(factor(sc_exprs_seurat_scaled@meta.data$zone)), 
       .f = function(obj_, name_) {
         obj_$df %>% 
           rename_at(.vars = vars(-Regulon),
                     .funs = function(x) {paste0(x, '_', name_)})
       }) %>%
  purrr::reduce(.f = left_join, by = 'Regulon')

#Collapse Data to avgNES and integrated p-value
stouffer_p <- function(p, weights) {metap::sumz(p = p, weights = weights)}
by_zone_aggregated_out <- by_zone_viper_df %>%
  melt(id.vars = 'Regulon') %>%
  tidyr::separate(col = 'variable', into = c('variable','zone'), sep = '_') %>%
  dcast(...~variable, value.var = 'value') %>%
  left_join(sc_exprs_seurat_scaled@meta.data %>%group_by(zone) %>% summarize(n=n()), 
            by = 'zone') %>%
  group_by(Regulon) %>%
  summarize('NES_avg' = weighted.mean(x = NES, w = sqrt(n)),
            'p.value' = stouffer_p(p = p.value, weights = sqrt(n))$p, #calling sumz directly with summarize cannot find variables...unclear why...
            'size' = mean(Size)) %>%
  mutate('fdr' = p.adjust(p.value, method = 'fdr'))
#by_zone_aggregated_out <- read_tsv('msviper_out.tsv')

#Look at the data
ggplot(data = by_zone_aggregated_out %>%
         mutate('identity' = ifelse(fdr<0.05, 'Significant','Not Significant')), 
       aes(x = NES_avg, y = -log10(fdr), color = identity)) +
  geom_point(alpha = 0.3) +
  ggrepel::geom_label_repel(data = . %>% filter(identity == 'Significant') %>%
                              arrange(fdr) %>% 
                              group_by(NES_avg > 0) %>%
                              dplyr::slice(1:3), 
                            aes(label = Regulon), show.legend = F, min.segment.length = unit(0, 'lines')) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_bw() +
  scale_color_manual(values = c('black','red') %>% set_names(c('Not Significant','Significant'))) +
  theme(legend.position = 'bottom',legend.title = element_blank()) +
  guides('color' = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  ylab('-log10(FDR)') + xlab('Average NES iRFP Bright/Dim')
