#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

echo $PWD
echo SPLIT : ${resources}

{

	java ${params.config.memory.splitStart} ${params.config.memory.splitMax} -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T SplitNCigarReads \
		-R ${resources}/human_g1k_v37.fasta \
		-I ${bam} \
		-o split.bam \
		${params.config.split.splitArgs}
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
