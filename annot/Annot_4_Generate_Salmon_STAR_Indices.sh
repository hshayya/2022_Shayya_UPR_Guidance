#Obtained sequences for fluorophores/transgenic sequences of interest
#gfp = EGFP from Clontech. https://www.addgene.org/vector-database/2485/
#tdtomato = from the plasmid that was used to make the Ai9 mice https://www.addgene.org/22799/
#LacZ from Ron Yu Addgene 32642 Ref. 2004 Neuron Paper (https://www.addgene.org/32642/sequences/). OMIT IRES AND TAU SEQUENCE TO AVOID homogeny to ires-tau-GFP etc.

#Create Combined fasta
cat IS_merged_CodingTranscriptsHS_transcripts.fasta egfp_tdtom_lacz.fa > IS_merged_CodingTranscriptsHS_transcripts_egfp_tdtom_lacz.fasta

#Create Salmon Index
salmon --no-version-check index -t IS_merged_CodingTranscriptsHS_transcripts_egfp_tdtom_lacz.fasta -i IS_merged_CodingTranscriptsHS_transcripts_egfp_tdtom_lacz_mm10_k31
 
#Create STAR Index
STAR   --runMode genomeGenerate   --runThreadN 20   --genomeDir /media/storageA/hani/alignment/STARDir_Coding75nt   --genomeFastaFiles /seq/mm10/iGenome/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa      --sjdbGTFfile IS_merged_CodingTranscriptsHS.gtf   --sjdbOverhang 74

#genome not included due to size. mm10 UCSC-style chrom numbers, from iGenomes.