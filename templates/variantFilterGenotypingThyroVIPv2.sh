#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	java ${params.config.memory.variantFilteringStart} ${params.config.memory.variantFilteringMax}  -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		-R ${resources}/human_g1k_v37.fasta \
		-V ThyroVIPv2.genotypes.raw.vcf \
		-o ThyroVIPv2.genotypes.filter.vcf \
		${params.config.variant_filtering.variantFilterGenotypingThyroVIPv2Args}
	ev=\$?

	sleep 60

} 2>&1

script_exit
