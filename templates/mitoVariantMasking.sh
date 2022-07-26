#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	java ${params.config.memory.variantFilteringStart} ${params.config.memory.variantFilteringMax} -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		-R ${resources}/human_g1k_v37.fasta \
		-V MTgenes.mt.filter.vcf \
		-o MTgenes.mt.filter.1000GHiConf.mask.vcf \
		--mask ${reference}/1000G_phase1.snps.high_confidence.b37.vcf \
		--maskName 1000GHiConf
	ev=\$?
	sleep 60
	

} 2>&1

script_exit
