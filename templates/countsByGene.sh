#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{
 
	$params.htseq_path/htseq-count \
		${params.config.countsByGene.htseqCountArgs} \
		${bam} \
		${resources}/Homo_sapiens.GRCh37.75.gtf > countsByGene.txt
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
