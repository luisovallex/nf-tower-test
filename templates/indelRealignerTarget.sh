#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	java ${params.config.memory.indelRealignerTargetStart} ${params.config.memory.indelRealignerTargetMax}  -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-R ${resources}/human_g1k_v37.fasta \
		-I ${bam} \
		-dt ${params.config.indel_realignment.downsamplingType} \
		-o target-intervals.list \
		-known ${resources}/Mills_and_1000G_gold_standard.indels.b37.vcf \
		-known ${resources}/1000G_phase1.indels.b37.vcf \
		${params.config.gatk_parallel.realignerTargetCreator}
	ev=\$?

	sleep 60
	

} 2>&1

script_exit
