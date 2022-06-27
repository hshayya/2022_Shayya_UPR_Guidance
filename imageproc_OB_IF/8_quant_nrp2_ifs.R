#Quantify Nrp2 Levels in M28 Ddit3 cKO tdtom vs. GFP glomeruli by IF

library(tidyverse)
library(magrittr)
library(reshape2)
setwd('/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/Analysis_Glomeruli')

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

#Newer Data appears to use #36ff00 whereas others uses #3cff00. Slightly Different Greens- probably some config updated somewhere. 
#We are NEVER directly comparing across animals so this shouldn't matter for any code anywhere.
list.files('/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/Analysis_Glomeruli/Extracted_Images/pxintz_outs/',
           full.names = T) %>%
  enframe() %>%
  mutate('imp' = str_extract(value, '([0-9]+\\-[0-9]+\\-[0-9]+_((B[0-9]+_[LR][1-5])|([0-9]+_4wkMor28))_(WT|Ctrl|Exp)[0-9]).*([LR](Lateral|Medial))'),
         'gfp_nm' = map_chr(value, .f = function(x) {colnames(pxintz_reader(x))[5]})) %>%
  group_by(imp, gfp_nm) %>%
  summarize(n=n())

#Parse Data and Convert to Log10
parsed_data <- 
  list.files('/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/Analysis_Glomeruli/Extracted_Images/pxintz_outs/',
             full.names = T) %>%
  map_dfr(.f = function(x) {
    df <- pxintz_reader(x)
    parsed_name <- str_split(x,'/')[[1]] %>% 
      .[length(.)] %>% str_split('\\.') %>% .[[1]] %>% .[1]
    df$parsed_name <- str_extract(parsed_name, '.*(WT|Ctrl|Exp)[0-9]+')
    df$loc <- str_extract(parsed_name, '[RL](Lateral|Medial)')
    df$section <- str_extract(parsed_name, '[0-9]+(?=_(gfp|tdt))')
    df$roi <- str_extract(parsed_name,'(gfp|tdt)$')
    
    df
  })

parsed_data <- parsed_data %>% 
  mutate_at(.vars = vars(magenta, red,green, `#0047ff`),
            .funs = log10)

#Red/Green fairly reproducible across animals & section/locations -> c/w genetically encoded
#Magenta staining more variable both across & within days
gfp_tdtom_vals <- c('gfp' = 'green','tdt' = 'red')
ggplot(data = parsed_data %>%
         melt(id.vars = c('X','Y','parsed_name','loc','section','roi'),
              measure.vars = c('magenta', 'red','green')), 
       aes(x = value, group = interaction(section, roi, variable), color = roi)) +
  geom_line(stat = 'density', alpha = 0.6) +
  facet_grid(rows = vars(parsed_name), cols = vars(variable), scales = 'free') +
  theme_bw() +
  scale_color_manual(values = gfp_tdtom_vals)

#"Quant Plots" which mimic the images generated
map(unique(parsed_data$parsed_name), .f = function(x) {
  df <- parsed_data %>% filter(parsed_name == x) %>%
    mutate('loc' = fct_relevel(loc, 'LLateral','LMedial', 'RLateral')) %>%
    group_by(loc) %>%
    mutate('facet_section' = fct_relevel(section, function(x) {as.character(sort(as.numeric(x)))}),
           'facet_section' = as.numeric(facet_section))
  
  ncol_ <- df %>% group_by(loc, section) %>% summarize() %>% summarize(n=n()) %>% pull(n) %>% max()
  
  ggplot(data = df,
         aes(x = magenta, color = roi)) +
    geom_line(stat = 'density') +
    theme_bw() +
    facet_grid(cols = vars(facet_section),
               rows = vars(loc)) +
    scale_color_manual(values = gfp_tdtom_vals) +
    ggtitle(x) +
    geom_text(data = . %>% group_by(loc, section, facet_section) %>% 
                summarize(),
              aes(label = section), x = -Inf, y = Inf, vjust = 1, hjust = 0, color = 'black',
              inherit.aes = F)
})

#Consider only sections with gfp and tdtom values
#Normalize Nrp2 signal in each region to the mean of the Nrp2 signal in GFP region
#The general trend is that there is less Nrp2 in the tdtomato glomeruli.
map(c('7-30-21_4250_4wkMor28_Exp1','9-16-21_4390_4wkMor28_Exp3'), .f = function(x) {
  ggplot(data = parsed_data %>%
           filter(parsed_name == x) %>%
           group_by(loc, section) %>%
           filter(length(unique(roi)) == 2) %>%
           mutate('corrected_magenta' = magenta - mean(magenta[roi == 'gfp'])) %>%
           ungroup() %>%
           mutate('loc' = fct_relevel(loc,'LLateral','LMedial','RLateral','RMedial'),
                  'section' = fct_relevel(section, function(x) {x[order(as.numeric(x))]})),
         aes(x = corrected_magenta, color = roi)) +
    geom_line(stat = 'density') +
    theme_bw() +
    facet_wrap(facets = vars(interaction(loc, section, lex.order = T))) +
    scale_color_manual(values = gfp_tdtom_vals) +
    ggtitle(x)
})

#Combined 4wk Data, for Supplements
#Again, we only consider images where there is BOTH a tdtomato and GFP roi
#We compute here the difference in log10 Nrp2 signal in tdt vs. GFP roi and plot as boxplot (combining n=2 biological replicates)
ggplot(data = parsed_data %>%
         filter(str_detect(parsed_name,'4wk')) %>%
         group_by(parsed_name, loc, section) %>%
         filter(length(unique(roi)) == 2) %>%
         summarize('corrected_magenta' = mean(magenta[roi == 'tdt']) - mean(magenta[roi == 'gfp'])) %>%
         ungroup(),
       aes(y = '1', x = corrected_magenta)) + #could use parsed name to show by animal!
  geom_point(position = position_jitter(height = 0.1)) +
  geom_boxplot(width = 0.2, fill = NA) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_bw() + ylab('') + 
  xlab('log10 Nrp2 signal in tdtom - GFP Glomeruli') +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

#Wilcox test for deviance from zero
parsed_data %>%
  filter(str_detect(parsed_name,'4wk')) %>%
  group_by(parsed_name, loc, section) %>%
  filter(length(unique(roi)) == 2) %>%
  summarize('corrected_magenta' = mean(magenta[roi == 'tdt']) - mean(magenta[roi == 'gfp'])) %>%
  ungroup() %>%
  pull(corrected_magenta) %>%
  wilcox.test()
