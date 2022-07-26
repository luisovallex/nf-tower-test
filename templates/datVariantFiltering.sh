#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

java ${params.config.memory.variantFilteringStart} ${params.config.memory.variantFilteringMax}  -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R ${resources}/human_g1k_v37.fasta \
	-V av2dat.genotypes.raw.vcf \
	-o av2dat.genotypes.filter.vcf \
	${params.config.variant_filtering.datVariantFilteringArgs}
ev=\$?

script_exit
