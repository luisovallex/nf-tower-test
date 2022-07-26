#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	# recalPre
	java ${params.config.memory.recalibrationStart} ${params.config.memory.recalibrationMax} -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-R ${resources}/human_g1k_v37.fasta \
		-I ${bam} \
		-o pre-recal-tbl.txt \
		-knownSites ${resources}/dbsnp_138.b37.vcf \
		-knownSites ${resources}/Mills_and_1000G_gold_standard.indels.b37.vcf \
		-knownSites ${resources}/1000G_phase1.indels.b37.vcf \
		${params.config.gatk_parallel.baseRecalibrator}
	ev=\$?

	

} 2>&1

script_exit
