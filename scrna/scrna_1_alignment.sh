#Build Cellranger Reference
#mm10 fasta not included for size reasons. 
cellranger mkref --genome mm10.ucsc_soria_pcdh --fasta mm10.fa --genes gene_soria_pcdh.gtf

#Cellranger v2.1 Alignment
cellranger count --fastqs=/path/to/fastqs --id=0718b --sample=SK003 --transcriptome=/path/to/reference

#Cellranger v2.2 Alignment
cellranger count --id=wt_rep2 --transcriptome=/path/to/reference --fastqs=/path/to/fastqs --sample=SR003 --localcores=16

#Cellranger v3.1 Alignment
cellranger count --id=WT --transcriptome=/path/to/reference --fastqs=/path/to/fastqs --sample=SJ003 --localcores=16