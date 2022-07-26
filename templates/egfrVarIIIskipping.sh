#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

echo "merged bam : ${bam}, metrics file : ${metrics}"
echo "genome ctat dir: ${resources}"

{
  # retrieve totalreads fromÂ  metricsfile
  totalreads="\$(cut -f 20 ${metrics} | tail -1)"
  echo "successfully retreived total reads"

  #run EGFR VIII : pyton file_name -i sampleName -n totalReads -f bamFilePath -a annotationFilePath eventName
  python /usr/local/setup/nextflow/bin/JR_METe14EGFRVIII.py -i $sample_name -n \${totalreads} -f ${bam} -a ${resources}/METe14EGFRVarIII_hg19_anno.csv -e EGFR_Var_III
  ev=\$?

  sleep 240
  
 } 2>&1

script_exit
