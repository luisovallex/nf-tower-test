#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	java ${params.config.memory.rnaqcStart} ${params.config.memory.rnaqcMax} -Djava.io.tmpdir=$params.javatmp -jar $params.rnaqc_path/RNA-SeQC_v1.1.8.jar \
		-n 1000 \
		-s "$sample_name|${bam}|NextSeq" \
		-t ${qcreference}/_RNASeQC_Refs/gencode.v7.annotation_goodContig.gtf \
		-r ${qcreference}/_RNASeQC_Refs/Homo_sapiens_assembly19.fasta \
		-o . \
		-strat gc \
		-gc ${qcreference}/_RNASeQC_Refs/gencode.v7.gc.txt
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
