#### R Script to Generate IS_merged_CodingTranscriptsHS.gtf ####

#Set Things Up
library(tidyverse)
library(magrittr)
setwd('/media/storageA/hani/alignment/mm10annotation_ensembl/')

####Import My HSannotations.gtf and partition into ORs from Ensembl and Ibarra Soria####
mygtf <- read.table('HSannotations.gtf',sep='\t',stringsAsFactors = F, header = F)

#Clean up the GTF to remove useless annotations. 
my_gtf_clean <- mygtf %>%
  mutate('gene'=str_extract(.$V9,'ENSMUSG.*?(?=;)'),
         'transcript' = str_extract(.$V9,'(?<=transcript_id ).*?(?=;)')) %>%
  filter(!is.na(transcript) & V3 %in% c('exon','CDS','5UTR','3UTR','start_codon','stop_codon'))

#Pull out list of Olfr Genes that have Ibarra Soria Transcripts
Olfr_genes <- my_gtf_clean %>%
  filter(str_detect(.$transcript,'CUFFORT')) %>%
  pull(gene)

#Pull out Ensembl Annotations for Transcripts in Olfr Genes
Olfr_gtf_ensembl<- my_gtf_clean %>%
  filter(gene %in% Olfr_genes & !str_detect(transcript, 'CUFFORT'))

#Pull out Ibarra Soria Transcripts in Olfr Genes
Olfr_gtf_IS <- my_gtf_clean %>%
  filter(str_detect(transcript, 'CUFFORT'))

#Write GTFs to file for analysis with Annot_2a_FindCDS_IS_transcripts.py
#write.table(Olfr_gtf_ensembl,"olfrs_ensembl.gtf",sep='\t',quote=F,row.names = F,col.names = F)
#write.table(Olfr_gtf_IS,"olfrs_IS.gtf",sep='\t',quote=F,row.names = F,col.names = F)

####Read in output from Annot_2a_FindCDS_IS_transcripts.py and Create Final GTFs####
#Read in my CDS annotations for the Ibarra_Soria Transcripts
new_annotations<-read.table('IS_ORseqs_CDSannotated_new.tsv',sep='\t',stringsAsFactors = F,header=F)
new_annotations %<>%
  mutate('gene'=str_extract(.$V9,'ENSMUSG.*?(?=;)'),
         'transcript' = str_extract(.$V9,'(?<=transcript_id ).*?(?=;)'))

####Add CDS annotations back to original GTF####
combined_annotations_clean<-rbind(my_gtf_clean,new_annotations)

#Write output with ONLY CODING Transcripts (Required for Metagene)####
#Get a list of transcripts that have a coding sequence
coding_transcripts <- combined_annotations_clean %>% 
  filter(V3 == 'CDS') %>%
  pull(transcript)

#Retain all annotations for coding transcripts 
gtf_coding <- combined_annotations_clean %>%
  filter(transcript %in% coding_transcripts) %>%
  select(-c(10:11))

#Write output
#write.table(gtf_coding,'IS_merged_CodingTranscriptsHS.gtf',sep='\t',quote=F,row.names = F,col.names = F)
