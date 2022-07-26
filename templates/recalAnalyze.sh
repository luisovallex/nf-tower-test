#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	# recalAnalyze
	java ${params.config.memory.recalibrationStart} ${params.config.memory.recalibrationMax}  -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T AnalyzeCovariates \
		-R ${resources}/human_g1k_v37.fasta \
		-before pre-recal-tbl.txt \
		-after post-recal-tbl.txt \
		-csv recalplots.csv
	ev=\$?

	sleep 60

} 2>&1

script_exit
