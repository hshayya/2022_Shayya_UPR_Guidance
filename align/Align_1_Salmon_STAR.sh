#Example Salmon Pseudo-Alignment. 
#libs_for_alignment.tsv: LibName\tR1_path\tR2_path (not included, can generate in R/python/excel)
#Index not included in dataset, follow code in ../annot to generate these
#Can use index +/- fluors. 
###### !!!!! PAPER USES +FLUORS for everything except OSN differentiation data. #####
#geneMap included on Zenodo, ./annot folder.
{
  while IFS=$'\t' read -ra line
  do
    echo Working on ${line[0]}
    salmon quant --no-version-check --libType A --seqBias --gcBias --index /path/to/Salmon/Index -1 ${line[1]} -2 ${line[2]} -o ./SalmonOut_${line[0]} -p 12 --geneMap tx_to_gene_name_gfp_tdtom_lacz.tsv --validateMappings
    echo ...done
  done
} < libs_for_alignment.tsv

#Example STAR Alignment
#Only did this for OmpCre Ddit3 libs
#libs_for_alignment.tsv: LibName\tR1_path\tR2_path (not included, can generate in R/python/excel)
#genomeDir not included in dataset, follow code in ../annot to make
{
  while IFS=$'\t' read -ra line
  do
    if [[ ${line[0]} =~ OmpCre ]]; then
      echo Working on File: ${line[0]}
      mkdir ./STARout_${line[0]}
      STAR --runMode alignReads --runThreadN 20 --genomeDir /path/to/genome/Dir --outFileNamePrefix ./STARout_${line[0]}"/"${line[0]}"_" --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMstrandField intronMotif --outFilterMultimapNmax 10 --readFilesIn ${line[1]} ${line[2]} --readFilesCommand zcat
    fi
  done
} < libs_for_alignment.tsv

#Example STAR post-processing
#First, find all bam files in any sub-directory of current working directory
rna_bams=()
while IFS=  read -r -d $'\0'; do
    rna_bams+=("$REPLY")
done < <(find . -iname *bam -type f -print0)

#Extract -q 30 reads and create bam index (all in the appropriate sub-dirs of current working directory)
for i in ${rna_bams[@]}; do
    echo Working on file $i
    my_path=$(echo $i | grep -Po ".*(?=_Aligned.sortedByCoord.out.bam)")
    samtools view -b -q 30 --threads 20 -o $my_path"_q30.bam" $i
    echo Indexing file $my_path"_q30.bam"
    samtools index $my_path"_q30.bam"
done