#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	java ${params.config.memory.variantCallingStart} ${params.config.memory.variantCallingMax}  -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T HaplotypeCaller  \
		-R ${resources}/human_g1k_v37.fasta \
		-I ${bam} \
		-o XAextension.2019Sep.genotypes.raw.vcf \
		-L ${resources}/XAextension.2019September.hg19.bed \
		${params.config.variant_calling.variantGenotypingArgs} \
		${params.config.gatk_parallel.haplotypeCaller}
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
