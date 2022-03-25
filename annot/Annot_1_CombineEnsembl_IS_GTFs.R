setwd("/media/storageA/hani/alignment/mm10annotation_ensembl/")

#GTF from Ensembl v92 pass GTF2.2 requirements. Fix this.

#Import GTF Files From Ensembl and Ibarra-Soria Plos Genetics Paper
ISframe<-read.table('journal.pgen.1004593.s017.GTF',sep='\t',stringsAsFactors = F,header=F,quote='') #ESSENTIAL to use quote='' otherwise it ignores quotations and messes everything up
ensemblframe<-read.table('Mus_musculus.GRCm38.92.gtf',sep='\t',stringsAsFactors = F,header = F,quote='')

#Adjust the V3 Column of Ensembl Frame to comply with specifications of GTF2
ensemblframe$V3[which(ensemblframe$V3 == "five_prime_utr")]<-"5UTR"
ensemblframe$V3[which(ensemblframe$V3 == "three_prime_utr")]<-"3UTR"

#All entries must have gene_id, transcript_id and gene_name in V9
#This is not the case in initial Ensembl GTF
length(which(!grepl("gene_id",ensemblframe$V9))) #everything has a gene_id
length(which(!grepl("transcript_id",ensemblframe$V9))) #53801 genes do not have transcript_id's (REQUIRED)
length(which(!grepl("gene_name",ensemblframe$V9))) #everything has a gene_name. Convenient but not required

#Add transcript_id "" to entries that are missing transcript ID
attributeslist<-strsplit(ensemblframe$V9,split=';')
correct_annotations<-function(annotation_vector){
  find_geneID<-which(grepl("gene_id",annotation_vector))[1]
  find_transcript<-which(grepl("transcript_id",annotation_vector))[1]
  find_genename<-which(grepl("gene_name",annotation_vector))[1]
  
  if(is.na(find_transcript)) {
    mystring<-paste(annotation_vector[find_geneID], " transcript_id \"\"", annotation_vector[find_genename], sep = ';')
    return(paste(mystring,";",sep=''))
  }
  else {
    mystring<-paste(annotation_vector[find_geneID],annotation_vector[find_transcript], annotation_vector[find_genename],sep=';')
    return(paste(mystring,';',sep=''))
  }
}

ensemblframe$V9<-unlist(lapply(attributeslist,correct_annotations))

#Ibarra-Soria uses "CUFFORG" gene_ids. Need to map back to Ensembl gene_ids using gene_names that they have in there
#Generate the gene_name -> gene_id lookup table
library(stringr)
annotatedgenes<-str_match(ISframe$V9,"gene_name (.*?);")[,1] #Its quicker to do in two steps than try to figure out a better regex
annotatedgenes<-unique(substr(annotatedgenes,12,nchar(annotatedgenes)-2)) #strip the quotation marks too

library(biomaRt) 
ensembl = useMart("ensembl", host = 'http://apr2018.archive.ensembl.org', dataset = 'mmusculus_gene_ensembl')

lookuptable<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                   filters = c("external_gene_name"),
                   values= annotatedgenes,
                   mart = ensembl)

annotatedgenes[which(!(annotatedgenes %in% lookuptable$external_gene_name))] #This works well, there are a few genes for which there is no Ensembl Annotation 
length(which(duplicated(lookuptable$external_gene_name))) #there are two entries where Ensembl assigned two gene names to one gene.
lookuptable<-lookuptable[!duplicated(lookuptable$external_gene_name),] #removes the second entry for each somewhat arbitrarily to avoid problems down the road. I looked and this doesn't seems super relevant- just Olfr1073-ps1 and Olfr290.

#Prep Ibarra-Soria GTF to be merged with Larger Ensembl GTF. 
#Keep only reconstructed Transcripts
toappend<-subset(ISframe,ISframe$V2=="reconstructed_transcript")

#Map Gene_names -> gene_ids
toappend$external_gene_name<-str_match(toappend$V9,"gene_name(.*?);")[,1]
toappend$external_gene_name<-substr(toappend$external_gene_name,12,nchar(toappend$external_gene_name)-2) #strip quotations
library(plyr)
toappend<-join(toappend,lookuptable,by="external_gene_name",type="left")
unique(toappend$external_gene_name[is.na(toappend$ensembl_gene_id)]) #12 genes with no ensembl annotation. 
toappend<-toappend[complete.cases(toappend),] #removes 103 rows with annotations belonging to those 12 genes. Would not be analyzing anyways. 

toappend$restofannotation<-str_extract(toappend$V9,"transcript_id.*")
toappend$V9<-paste("gene_id \"",toappend$ensembl_gene_id, "\"; ",toappend$restofannotation,sep='')
appendedframe<-rbind(ensemblframe,toappend[,c(1:9)])

#Lastly, fix chromosomes
chromstokeep<-as.character(c(1:19,"X","Y","MT"))
appendedframe<-appendedframe[(appendedframe$V1 %in% chromstokeep),] #we don't have the unlocalized contigs in our genome fasta! no need to keep.
appendedframe$V1<-paste0("chr",appendedframe$V1)
appendedframe$V1[appendedframe$V1=="chrMT"]<-"chrM"

#write.table(appendedframe,"HSannotations.gtf",sep='\t',quote=F,row.names = F,col.names = F) 

#GTF seems to now mostly pass the GTF2 validation script. 