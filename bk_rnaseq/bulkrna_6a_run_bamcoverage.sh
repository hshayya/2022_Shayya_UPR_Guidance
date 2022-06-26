#Example script showing how bamCoverage was run to generate normalized bw files for OmpCre Ddit3 Experiment
#See bulkrna_6_Viz_OmpCreDdit3_CoveragePlots.R for generation of the parameter file used here
#Tab sep contents: bam file path | output file path prefix | 1/sizefactor 
{
  while IFS=$'\t' read -ra line
  do
    echo Working on: ${line[0]}
    echo Forward Strand...
    bamCoverage -b ${line[0]} -o ${line[1]}"fwd.bw" --filterRNAstrand forward --scaleFactor ${line[2]} -p 10 --binSize 1
    echo Reverse Strand...
    bamCoverage -b ${line[0]} -o ${line[1]}"rev.bw" --filterRNAstrand reverse --scaleFactor ${line[2]} -p 10 --binSize 1
  done
} < /path/to/parameter/file.tsv