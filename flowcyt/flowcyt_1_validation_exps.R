library(tidyverse)
library(reshape2)
library(patchwork)

#Set Everything Up

setwd('./')
csv_path <- './flowcyt'

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

#Load Validation Data
validation_data <- 
  list.files(csv_path,
             full.names = T) %>%
  .[str_detect(.,'12_Feb_21')] %>%
  map_dfr(.f = function(x) {
    fname <- str_split(x,'/')[[1]] %>% .[length(.)]
    read_csv(x) %>%
      dplyr::select('GFP' = 'FL17-Area', 'iRFP' = 'FL37-Area') %>%
      mutate('population' = str_extract(fname,'iRFP\\+(GFP(Bright|Dim))?(?=.csv$)'),
             'date' = str_extract(fname, '[0-9]{2}_[A-Za-z]{3}_(20|21)'),
             'tube' = paste(str_extract(fname, '[0-9]{3}(?=_Exp)'), 
                            str_extract(fname,'(?<=GFP\\-iRFP ).*?(?=_12)'),
                            sep='_'))
  })

#Compute Means/SDs for iRFP+ population in each tube. 
mean_sd_table <- validation_data %>%
  filter(population == 'iRFP+') %>%
  melt(id.vars = 'tube',
       measure.vars = c('GFP','iRFP')) %>%
  group_by(tube, variable) %>% 
  summarize(mean = mean(value), 
            sd = sd(value))

#Correct the data by "centering" (subtract iRFP+ mean)
validation_data_out <- validation_data %>% 
  mutate('cell_id' = seq(1, length(tube))) %>%
  melt(measure.vars = c('GFP','iRFP')) %>%
  left_join(mean_sd_table, by = c('tube','variable')) %>%
  mutate('centered_value' = value-mean) %>%
  dplyr::rename('raw_value' = 'value') %>%
  melt(id.vars = c('population','tube','cell_id','variable'),
       measure.vars = c('raw_value', 'centered_value'),
       variable.name = 'analysis_type') %>%
  dcast(...~variable, value.var = 'value')

#Plot raw and corrected data, showing voltage shifts resolved.
ggplot(data = validation_data_out %>%
         group_by(tube) %>%
         filter(cell_id %in% sample(c(min(cell_id):max(cell_id)), 1000)) %>%
         ungroup() %>%
         mutate('tube' = fct_recode(tube, '0V_0V' = '159_0V',
                                    '+50V_+50V' = '160_+50V',
                                    '-50V_-50V' = '161_-50V',
                                    '-50V_+50V' = '162_-50V +50V',
                                    '+50V_-50V' = '163_-50V +50V'),
                'analysis_type' = fct_recode(analysis_type, 
                                             'Raw' = 'raw_value',
                                             'Corrected' = 'centered_value')),
       aes(x = iRFP, y = GFP, color = tube)) +
  geom_point(alpha = 0.1) +
  facet_wrap(facets = vars(analysis_type), scales = 'free') +
  theme_bw() +
  guides(color = guide_legend(nrow = 2, override.aes = list(alpha = 1))) +
  theme(legend.position = 'bottom', legend.title = element_blank())

#Load OR-iGFP iRFP measurement data
load_data <- list.files(csv_path,
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

#Show that raw iRFP and GFP values for OR-iGFP experiments are contained within the raw ranges tested by the PMT experiment
#This suggests that our PMT experiment will provide a good estimate of technical noise after the correction approach
ggplot(data = validation_data_out %>%
         filter(analysis_type == 'raw_value' & population %in% c('iRFP+','iRFP+GFPBright')) %>%
         rename(date = tube) %>%
         list(load_data) %>% 
         set_names(c('pmt_experiments','or_experiments')) %>%
         imap(.f = function(.x,.y) {
           melt(.x, 
                id.vars = c('population','date'), 
                measure.vars = c('GFP','iRFP')) %>%
             group_by(population, date, variable) %>%
             summarize(mean = mean(value),
                       min = mean - sd(value),
                       max = mean + sd(value)) %>%
             ungroup() %>%
             mutate('experiment' = rep(.y, nrow(.)))
         }) %>% 
         bind_rows() %>%
         melt(measure.vars = c('mean','min','max'),
              variable.name = 'measurement_type') %>%
         dcast(...~variable+measurement_type, value.var = 'value') %>%
         mutate(experiment = fct_recode(experiment, 'OR-iGFP Experiments' = 'or_experiments',
                                        'Technical Controls' = 'pmt_experiments'),
                experiment = fct_relevel(experiment,'Technical Controls')),
       aes(x = iRFP_mean, y = GFP_mean, color = experiment)) +
  geom_point() +
  geom_errorbar(aes(ymin = GFP_min, ymax = GFP_max), alpha = 0.3) +
  geom_errorbarh(aes(xmin = iRFP_min, xmax = iRFP_max), alpha = 0.3) +
  theme_bw() + xlab('iRFP') + ylab('GFP') +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2))
