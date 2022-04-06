library(tidyverse)
library(reshape2)
library(readxl)
library(patchwork)

#Example showing approach to blinding OE IF images for cell counts in M71/Mor28 Perk Mice
#These experiments were done in batches, with each batch having a mix of WT/Het/cKO mice
#Downsampled images used for all counts -> faster to load/can still easily see +ve vs. -ve cells

gen_letters_indefinitely <- function(len_, case = 'upper') {
  letter_ <- switch(case,'upper' = LETTERS, 'lower' = letters)
  iterations <- ceiling(len_/length(letter_))
  out <- unlist(map(.x = seq(1,iterations),
                    .f = function(x) {do.call(paste0, map(seq(1,x), function(a) {letter_}))}))
  out[1:len_]
}

#Map images -> random blinded file names
set.seed(2)
m28_p5_stain_lookuptable_11_4 <- 
  list.files('path/to/downsampled/images',
             pattern = '.*tif', full.names = T) %>%
  list() %>%
  set_names('path') %>% 
  as.data.frame() %>%
  mutate('filename' = map_chr(path, .f = function(x) {
    str_split(x,'/')[[1]] %>% .[length(.)]
  })) %>%
  tidyr::extract(col = filename, into = c('Animal','Genotype','Slide','Section'),
                 regex = '([0-9]+-[0-9]+-[0-9]+_B[0-9]+_[L|R][0-9](?:[L|R][0-9])?).*(Ctrl|Exp|WT)(?:[0-9]+)?_([0-9]+)_([A-D])',
                 remove = F) %>% #Files always named animal_genotype_slide_section. Animal = Date_Toe([LR][1-5]). Gt WT = WT, Het = Ctrl, cKO = Exp
  mutate('random_code' = sample(gen_letters_indefinitely(len_ = nrow(.),
                                                         case = 'lower'),
                                replace = F),
         'new_file' = paste0('/path/to/blinded/analysis/location',
                             random_code, '.tif'))

#Save the Lookuptable
#write_tsv(m28_p5_stain_lookuptable_11_4, 
#          path = '/path/to/lookup/table')

#Copy the Images with the blinded names
walk2(.x = m28_p5_stain_lookuptable_11_4$path,
      .y = m28_p5_stain_lookuptable_11_4$new_file,
      .f = function(old_, new_) {
        print(paste('Copying', old_, 'to', new_, sep = ' '))
        file.copy(from = as.character(old_),
                  to = as.character(new_), overwrite = F)
      })

#Read into R for plotting/additional analysis
m28_p5_replicate1 <- m28_p5_stain_lookuptable_11_4 %>%
  left_join(read_excel('path/to/blinded/counts/excel/file'), 
            by = c('random_code' = 'Image')) %>%
  mutate('ratio' = `tdtom+M28+`/`GFP+M28+`) 
#excel file with columns Image = random_code, tdtom+M28-,	tdtom+M28+,	GFP+M28-,	GFP+M28+,	GFP+tdtom+M28+
#one row for each OE section