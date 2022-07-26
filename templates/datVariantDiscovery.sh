#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	java ${params.config.memory.variantCallingStart} ${params.config.memory.variantCallingMax}  -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T HaplotypeCaller  \
		-R ${resources}/human_g1k_v37.fasta \
		-I ${bam} \
		-o av2dat.genotypes.raw.vcf \
		-L ${resources}/AV2DATCap.gte5UniqSmpl_1000GHiConf_CADDgt10.v7M.bed \
		${params.config.variant_calling.datVariantDiscoveryArgs} \
		${params.config.gatk_parallel.haplotypeCaller}
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
