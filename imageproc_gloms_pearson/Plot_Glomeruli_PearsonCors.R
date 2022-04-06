library(tidyverse)
library(magrittr)
library(reshape2)
setwd('/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Experiments/Analysis_Glomeruli/Zstacks_Confocal_Quant/')

#Example code to read .pxintz files measuring red/green pixel intensities in 20x glomerular images -> compute per-glomerulus pearson correlations -> plot.

#'Parse .pxintz file
#'@param fpath (str) specifying disk location of .pxintz file to parse
#'@output (tibble) with parsed values
pxintz_reader <- function(fpath) {
  con <- file(fpath, open = 'r')
  
  output <- list() #hold the pixel values
  current_title <- NA
  current_z <- NA
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (substr(oneLine,1,1) == '#') {
      current_z <- substr(oneLine, 2, nchar(oneLine))
      output[[current_z]] <- list()
    }
    else if (substr(oneLine,1,1) == '>') {
      current_title <- substr(oneLine, 2, nchar(oneLine))
    }
    else {
      output[[current_z]][[current_title]] <- as.numeric(strsplit(oneLine, ',')[[1]])
    }
  } 
  close(con)
  
  #Map through and add the names
  imap_dfr(output, function(ls_, name_) {
    out <- as_tibble(ls_) 
    out$Z <- rep(as.numeric(name_), nrow(out))
    return(out)
  })
}

#Load Files and Parse
quant_files <- list.files('pxintz_outs/',
                          full.names = T)

pb <- progress_estimated(length(quant_files))

all_quants <- quant_files %>%
  map_dfr(.f = function(x) {
    df <- pxintz_reader(x)
    name_ <- x %>% str_split('/') %>% .[[1]] %>% .[length(.)] %>% 
      str_extract('.*(?=\\.pxintz)')
    
    out <- df %>% mutate('name' = rep(name_, nrow(.)))
    
    pb$tick()$print()
    return(out)
  }) 

all_quants <- all_quants %>% 
  mutate('#3cff00' = ifelse(is.na(`#3cff00`), `#36ff00`, `#3cff00`)) 
#different GFP channel names on some images... this fixes it
#note that we NEVER make quantiative comparisons BETWEEN glomeruli, its always correlations red/green WITHIN glomeruli
#differences in gfp channel parameters between images are thus irrelevent.

#Convert to Log10
by_OB_by_glom <- all_quants %>%
  group_by(name) %>%
  summarize('cor' = cor.test(log10(red), log10(`#3cff00`), method =  'pearson')$estimate) %>%
  ungroup() 

#Viz the correlations by Glomerulus (Group ~ Genotype)
ggplot(data = by_OB_by_glom %>%
         tidyr::extract(col = 'name', into = c('Animal','Genotype','Glomerulus'),
                        regex = '([0-9]+\\-[0-9]+(?:\\-[0-9]+)?_B[0-9]+_[LR1-9]+).*(WT|Ctrl|Exp).*((?:L|R)(?:Lateral|Medial|medial|meidla))') %>%
         arrange(Genotype, Animal, Glomerulus) %>%
         mutate('Glomerulus' = fct_collapse(Glomerulus, 'LMedial' = c('LMedial','Lmedial'),
                                            'RMedial' = c('RMedial','Rmedial','Rmeidla'))) %>%
         mutate('Genotype' = fct_recode(Genotype,'Het' = 'Ctrl','cKO' = 'Exp'),
                'Genotype' = fct_relevel(Genotype,'WT','Het')),
       aes(x = Genotype, color = Genotype, y = cor)) +
  #  geom_point(aes(group = Animal),
  #             position = position_jitterdodge(jitter.width = 0.05,
  #                                             dodge.width = 0.5),
  #             alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.15),
             alpha = 0.5) +
  geom_boxplot(width = 0.3, fill = NA, outlier.shape = NA) +
  geom_text(data = . %>% 
              lm(formula = cor~Genotype) %>%
              anova() %>%
              .['Genotype','Pr(>F)'] %>%
              {data.frame('p' = .)}%>%
              melt(),
            aes(label = paste0(variable, '=',formatC(value, digits = 2, format = 'e'))), 
            x = Inf, y=Inf, hjust=1, vjust =1, inherit.aes = F) +
  geom_text(data = . %>% group_by(Genotype) %>% 
              summarize(n=n()) %>%
              melt(id.vars = 'Genotype'), 
            aes(label = paste0(variable, '=', value),
                y = 0.8), show.legend = F) +
  theme_bw() + expand_limits(y = 0.9) +
  xlab('') + ylab('Pearson Correlation Red/Green') + ggtitle('p5 Mor28 Ddit3 Glomeruli')
