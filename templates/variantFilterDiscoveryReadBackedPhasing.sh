#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

echo "XXX check"

{

	java ${params.config.memory.variantFilteringStart} ${params.config.memory.variantFilteringMax} -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T ReadBackedPhasing \
		-R ${resources}/human_g1k_v37.fasta \
		-I ${bam} \
		-V IllmRNAAccessCap.discover.filter.vcf \
		-o IllmRNAAccessCap.discover.filter.phased.vcf \
		-enableMergeToMNP \
		--phaseQualityThresh 20.0 \
		--maxGenomicDistanceForMNP 10
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
