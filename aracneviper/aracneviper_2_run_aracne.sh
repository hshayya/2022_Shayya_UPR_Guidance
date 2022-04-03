#!/bin/bash

#Compute MI Threshold
java -Xmx40G -jar /media/storageA/hani/software/ARACNe-AP-master/dist/aracne.jar -e 272samples_OE_expression_matrix_forArachne.tsv -o Arachne_Out --tfs tf_list.txt --pvalue 1E-8 --seed 1 --calculateThreshold

#Run 100 bootstraps
for i in {1..100}; do
  echo Bootstrapping Iteration $i of 100...
  java -Xmx40G -jar /media/storageA/hani/software/ARACNe-AP-master/dist/aracne.jar -e 272samples_OE_expression_matrix_forArachne.tsv  -o Arachne_Out --tfs tf_list.txt --pvalue 1E-8 --seed $i --threads 12
done

#Consolidate Bootsraps
java -Xmx40G -jar /media/storageA/hani/software/ARACNe-AP-master/dist/aracne.jar -o Arachne_Out --consolidate