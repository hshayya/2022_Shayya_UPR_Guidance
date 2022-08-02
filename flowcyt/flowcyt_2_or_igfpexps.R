library(tidyverse)
library(reshape2)
library(patchwork)

#Set things up
setwd('./')
csv_dir <- './'

#'Autoscale a Vector
#'@description Linearly scale x to the range [0,1] 
#'@param x (double) numeric vector to be scaled
#'@param saturated (double) percentage of events to consider saturated (set to Inf)
#'@param invert (bool) should input be inverted prior to scaling
#'@return (double) scaled and (possibly) inverted 
autoscale_LUT <- function(x, saturated = 0.35, invert = F) {
  if (invert) {
    x <- -x
  }
  out <- x-min(x)
  out <- out/quantile(out, 1-(saturated/100))
  out[out>1]<-Inf
  out
} 

#Load all the OR-iGFP Experiments
load_data <- list.files(csv_dir,
                        full.names = T) %>% 
  .[str_detect(.,'09_Nov_20|11_Nov_20|02_Dec_20|27_Jan_21|17_Feb_21|23_Feb_21|10_Mar_21|16_Mar_21|14_Sep_20|14_Oct_20|16_Jun_21')] %>% 
  .[!str_detect(., 'Neg|[Cc]ontrol')] %>% 
  .[!str_detect(.,'M23-ATF-1 GFP-iRFP_10_Mar_21')] %>%
  map_dfr(.f = function(x) {
    fname <- str_split(x,'/')[[1]] %>% .[length(.)]
    read_csv(x) %>%
      dplyr::select('GFP' = 'FL17-Area', 'iRFP' = 'FL37-Area') %>%
      mutate('population' = str_extract(fname,'iRFP\\+(GFP\\+)?(?=.csv$)'),
             'date' = str_extract(fname, '[0-9]{2}_[A-Za-z]{3}_(20|21)'),
             'tube' = str_extract(fname,'.*(?=iRFP\\+(GFP\\+)?.csv$)'))
  })

#The below code computes mean/SD for GFP and iRFP EACH DAY. 
#So while we record tube #'s for each day, the underlying assumption is that PMT/technical variation is happening across days more so than between tubes
#The EXCEPTION is on 3/16/21, where we changed the PMT settings for the second tube (YFP was a little brighter than expected initially...)
#Easiest way to make this work is to split day from 3/16/21 into 2 "days" so you can just run the pipleine normally after this
load_data <- load_data %>%
  mutate('date' = ifelse(date == '16_Mar_21', paste(date,as.character(as.numeric(factor(tube))), sep='_'), date)) 

#Compute mean/SD for GFP and iRFP measurements in iRFP+ population EACH DAY.
#See above for why we do this per day not per tube.
sd_table <- load_data %>% 
  filter(population == 'iRFP+') %>%
  melt(id.vars = 'date',
       measure.vars = c('GFP','iRFP')) %>%
  group_by(date, variable) %>%
  summarize('mean' = mean(value),
            'sd' = sd(value))

#Normalize all data by the iRFP+ population for that day. 
#Centering approach, so just subtract mean!
output_data <- load_data %>% 
  group_by(population,date,tube) %>%
  mutate('cell_id' = seq(1, length(population))) %>% #careful with this. by TUBE!
  melt(measure.vars = c('GFP','iRFP')) %>%
  left_join(sd_table, by = c('date','variable')) %>%
  mutate('norm_value' = value-mean) %>%
  dcast(cell_id + date + tube + population ~ variable, value.var = 'norm_value')

#OR/date lookups
date_lookup <- bind_rows(list(
  data.frame('date' = c('11_Nov_20','02_Dec_20'), 'OR' = rep('P2',2)),
  data.frame('date' = c('27_Jan_21', '23_Feb_21'), 'OR' = rep('Mor28',2)),
  data.frame('date' = c('09_Nov_20', '14_Sep_20', '14_Oct_20'), 'OR' = rep('M71',3)),
  data.frame('date' = c('17_Feb_21', '10_Mar_21'), 'OR' = rep('Mor23',2)),
  {sd_table %>% 
      filter(str_detect(date, '16_Mar_21|16_Jun_21')) %>% 
      group_by(date) %>% 
      summarize() %>%
      mutate('OR' = rep('FishOR',nrow(.)))
  })
)

#Prepare the Data for Plotting (only need iRFP+GFP+, add OR mapping, compute "differentiation axis")
output_data_toplot <- output_data %>%
  filter(population == 'iRFP+GFP+') %>%
  left_join(date_lookup, by = 'date') %>%
  group_by(OR) %>%
  nest() %>%
  mutate('data' = map2(OR, data, function(k,v) {
    if (k == 'FishOR') {mutate(v, 'differentiation_axis' = autoscale_LUT(GFP, invert = T))}
    else {mutate(v, 'differentiation_axis' = autoscale_LUT(GFP, invert = F))}
  })) %>%
  unnest(cols = 'data') %>%
  ungroup() %>%
  mutate('OR' = fct_relevel(OR,'FishOR','M71','Mor23','P2','Mor28'))

#Additional Preparations for Plotting
#   (a) Subsample 400 random cells for each OR 
#   (b) Add a "merged" ClassI/M71/Mor23 panel by duplicating relevant data under a second facet variable
output_data_toplot_formerged <- 
  output_data_toplot %>% group_by(OR) %>%
  slice_sample(n = 400) %>%
  ungroup() %>%
  left_join(data.frame('OR' = c('FishOR','Mor23','M71','P2','Mor28'),
                       'zone' = c(rep('Zone1',3), 'Zone2','Zone5')),
            by = 'OR') %>%
  mutate('facet_var' = OR)

output_data_toplot_formerged <- 
  bind_rows(output_data_toplot_formerged,
            output_data_toplot_formerged %>% filter(zone == 'Zone1') %>% mutate('facet_var' = rep('Dorsal_Merged'))) %>%
  mutate('facet_var' = fct_relevel(facet_var, 'FishOR','M71','Mor23','Dorsal_Merged','P2'))

#Point Plot
new_color_scale_2_17 <- scales::hue_pal()(7)[c(1,3,5,6,7)] %>% set_names(c('M71','Class I','Mor23','P2','Mor28'))
ggplot(data = output_data_toplot_formerged %>%
         mutate('OR' = fct_recode(OR, 'Class I' = 'FishOR'),
                'facet_var' = fct_recode(facet_var, 'Class I' = 'FishOR')) %>%
         mutate('facet_var' = fct_relevel(facet_var, 'Class I','M71','Mor23','Dorsal_Merged')),
       aes(x = iRFP, y = differentiation_axis, color = OR)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(facets = vars(facet_var), nrow = 2) + #SL prefers this way, I like ~OR but either is fine
  ylab('Differentiation Axis') + xlab('Corrected iRFP') +
  xlim(c(-150,150)) + guides(color = guide_legend(override.aes = list(alpha=1))) +
  scale_color_manual(values = new_color_scale_2_17)

#Loess Smoothed iRFP ~ Differentiation Axis Plot
#Note that this plot uses ALL cells (contrast w/ point plot above where subsample 400 cells/OR)
ggplot(data = output_data_toplot %>%
         mutate('OR' = fct_recode(OR, 'Class I' = 'FishOR')),
       aes(x = differentiation_axis, y = iRFP, color = OR, fill = OR, group = OR)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_smooth(method = 'loess', formula = y~x, se = T, alpha = 0.3) +
  theme_bw() +
  xlab('Differentiation Axis') + ylab('Corrected iRFP') +
  scale_color_manual(values = new_color_scale_2_17) +
  scale_fill_manual(values = new_color_scale_2_17)

#QC for Supplementals: Show there isn't a counfound by day
ggplot(data = output_data_toplot %>%
         group_by(OR, date) %>% 
         dplyr::slice_sample(n = 200) %>%
         group_by(OR) %>%
         mutate('Sample' = as.character(as.numeric(factor(date))),
                OR = fct_recode(OR, 'Class I' = 'FishOR'),
                OR = fct_relevel(OR, 'Class I', 'M71','Mor23','P2','Mor28')),
       aes(x = iRFP, y = differentiation_axis, color = Sample)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_bw() +
  facet_wrap(facets = vars(OR), nrow = 2) + ylab('Differentiation Axis')

#Load in the Perk Hemizygous Data
all_perkhemi_data <- 
  list.files(csv_dir,
             full.names = T) %>%
  .[str_detect(., '31_Mar_21|16_Aug_21|20_Sep_21') & !str_detect(.,'Gating|M28')] %>%
  map(.f = function(x) {
    fname <- str_split(x,'/')[[1]] %>% .[length(.)]
    read_csv(x) %>%
      dplyr::select('GFP' = 'FL17-Area', 'iRFP' = 'FL37-Area') %>%
      mutate('population' = str_extract(fname,'iRFP\\+(GFP\\+)?(?=.csv$)'),
             'date' = str_extract(fname, '[0-9]{2}_[A-Za-z]{3}_(20|21)'),
             'tube' = str_extract(fname,'.*(?=iRFP\\+(GFP\\+)?.csv$)'))
  }) %>%
  .[map_lgl(., .f = function(x) {nrow(x) > 1})] %>% #remove the populations without any cells 
  bind_rows() %>%
  mutate('sample_id' = ifelse(str_detect(tube,'Control|CTRL'),
                              'Perk_WT','Perk_Het')) %>%
  group_by(date) %>%
  mutate('iRFP' = iRFP - mean(iRFP[sample_id == 'Perk_WT']), #compute Het-WT iRFP within replicate
         'facet' = 'OmpGFP',
         'OR' = 'OmpGFP',
         'genotype' = sample_id,
         'differentiation_axis' = autoscale_LUT(GFP)) %>% #autoscaling within replicate. Reasonable since not correcting GFP within replicate as we did for OR-iGFP (there autoscale COMBINED data)
  group_by(date, genotype) %>% 
  dplyr::slice_sample(n = min(group_size(.))) %>% #equal weights for each replicate
  ungroup()

#Plot Perk Het vs. WT iRFP Levels on top of OR-iGFP Exp Loess Results
ggplot(data = output_data_toplot %>%
         dplyr::select(OR, differentiation_axis, iRFP) %>%
         mutate('facet' = rep('OR',nrow(.)),
                'OR' = fct_recode(OR, 'Class I' = 'FishOR'),
                'sample_id' = OR,
                'genotype' = rep('Perk_WT', nrow(.))),
       aes(x = differentiation_axis, y = iRFP, group = sample_id, 
           color = OR, fill = OR)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_smooth(method= 'loess', geom = 'line', alpha = 0.2, size  =1) +
  geom_smooth(method = 'loess', formula = y~x, se = T, alpha = 0.1, linetype = 0) +
  geom_smooth(data = dplyr::select(dplyr::slice_sample(all_perkhemi_data, 
                                                       n=nrow(output_data_toplot)), #subsample Omp-GFP data, lots of cells otherwise!
                                   OR,differentiation_axis, iRFP,facet, sample_id, genotype) %>%
                mutate('genotype' = fct_relevel(genotype,'Perk_WT')),
              aes(linetype = genotype),
              method = 'loess',
              alpha = 0.3) +
  theme_bw() +
  xlab('Differentiation Axis') + ylab('Corrected iRFP') +
  scale_color_manual(values = c(new_color_scale_2_17,c('OmpGFP' = 'black')),
                     breaks = names(c(new_color_scale_2_17,c('OmpGFP' = 'black')))) +
  scale_fill_manual(values = c(new_color_scale_2_17,c('OmpGFP' = 'black')),
                    breaks = names(c(new_color_scale_2_17,c('OmpGFP' = 'black')))) +
  guides(color = guide_legend(title = 'Population'), fill = guide_legend(title = 'Population'),
         linetype = guide_legend(title = 'Genotype', override.aes = list(color = 'black', fill = NA))) +
  theme_bw() + xlab('Differentiation Axis') 

#Load in the Hspa Heterozygous Data
omp_hsp_atf5 <-
  list.files(csv_dir,
             full.names = T) %>%
  .[str_detect(., '27_Jul_22')] %>%
  map(.f = function(x) {
    fname <- str_split(x,'/')[[1]] %>% .[length(.)]
    read_csv(x) %>%
      dplyr::select('GFP' = 'FL17-Area', 'iRFP' = 'FL37-Area') %>%
      mutate('population' = str_extract(fname,'iRFP\\+(GFP\\+)?(?=.csv$)'),
             'date' = str_extract(fname, '[0-9]{2}_[A-Za-z]{3}_(20|21|22)'),
             'tube' = str_extract(fname,'.*(?=iRFP\\+(GFP\\+)?.csv$)'))
  }) %>%
  .[map_lgl(., .f = function(x) {nrow(x) > 1})] %>% #remove the populations without any cells 
  bind_rows() %>%
  mutate('sample_id' = str_extract(tube,'WT|Het'),
         'sample_id' = ifelse(sample_id == 'WT','Hspa5_WT','Hspa5_Het'),
         'sample_id' = fct_relevel(sample_id, 'Hspa5_WT'))

#Normalize to Hsp WT and compute differentiation axis
omp_hsp_toplot <- omp_hsp_atf5 %>%
  mutate('iRFP' = iRFP - mean(iRFP[sample_id == 'Hspa5_WT']),
         'facet' = 'OmpGFP',
         'OR' = 'OmpGFP',
         'genotype' = sample_id,
         'differentiation_axis' = autoscale_LUT(GFP))

set.seed(1)
ggplot(data = output_data_toplot %>%
         dplyr::select(OR, differentiation_axis, iRFP) %>%
         mutate('facet' = rep('OR',nrow(.)),
                'OR' = fct_recode(OR, 'Class I' = 'FishOR'),
                'sample_id' = OR,
                'genotype' = rep('Perk_WT', nrow(.))),
       aes(x = differentiation_axis, y = iRFP, group = sample_id, 
           color = OR, fill = OR)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_smooth(method= 'loess', geom = 'line', alpha = 0.2, size  =1) +
  geom_smooth(method = 'loess', formula = y~x, se = T, alpha = 0.1, linetype = 0) +
  geom_smooth(data = dplyr::select(dplyr::slice_sample(omp_hsp_toplot, 
                                                       n=nrow(output_data_toplot)), #subsample Omp-GFP data, lots of cells otherwise!
                                   OR,differentiation_axis, iRFP,facet, sample_id, genotype),
              aes(linetype = genotype),
              method = 'loess',
              alpha = 0.3) +
  theme_bw() +
  xlab('Differentiation Axis') + ylab('Corrected iRFP') +
  scale_color_manual(values = c(new_color_scale_2_17,c('OmpGFP' = 'black')),
                     breaks = names(c(new_color_scale_2_17,c('OmpGFP' = 'black')))) +
  scale_fill_manual(values = c(new_color_scale_2_17,c('OmpGFP' = 'black')),
                    breaks = names(c(new_color_scale_2_17,c('OmpGFP' = 'black')))) +
  guides(color = guide_legend(title = 'Population', order = 1), 
         fill = guide_legend(title = 'Population', order = 1),
         linetype = guide_legend(title = 'Genotype', override.aes = list(color = 'black', fill = NA))) +
  theme_bw() + xlab('Differentiation Axis') 
