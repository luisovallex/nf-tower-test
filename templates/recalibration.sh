#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	# recalibration
	java ${params.config.memory.recalibrationStart} ${params.config.memory.recalibrationMax} -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T PrintReads \
		-R ${resources}/human_g1k_v37.fasta \
		-I ${bam} \
		-BQSR pre-recal-tbl.txt \
		-o calibrated.bam \
		${params.config.gatk_parallel.printReads}
	ev=\$?

	sleep 60

	

} 2>&1

script_exit
