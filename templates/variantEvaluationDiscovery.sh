#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	java ${params.config.memory.variantEvaluationStart} ${params.config.memory.variantEvaluationMax}  -Djava.io.tmpdir=$params.javatmp -jar $params.gatk_path/GenomeAnalysisTK.jar \
		-T VariantEval \
		-R ${resources}/human_g1k_v37.fasta \
		--dbsnp ${resources}/dbsnp_138.b37.vcf \
		-eval IllmRNAAccessCap.discover.filter.vcf  \
		-o IllmRNAAccessCap.discover.filter.eval.txt \
		${params.config.gatk_parallel.variantEval}
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
