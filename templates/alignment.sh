#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	echo "runid = ${runid}"
	echo "sample_name = ${sample_name}"
	echo "sample_id = ${sample_id}"
	echo "lane_id = ${lane_id}"
	echo "readFilesIn = ${read1} ${read2}"
	echo "uuid = ${params.uuid}"
	echo "reference = ${reference}"
	echo "params_reference = ${params.reference}"

	echo "starting STAR at `date`"

	echo "making temporary directory ${params.javatmp}/${params.uuid}"
	mkdir -p ${params.javatmp}/${params.uuid}

	echo "removing temp dir in case of a retry..."
	rm -rf ${params.javatmp}/${params.uuid}/STAR-${runid}-${sample_name}.${lane_id}

	$params.star_path/STAR \
		--genomeLoad ${params.config.alignment.genomeLoad} \
		--genomeDir ${params.reference} \
		--readFilesIn ${read1} ${read2} \
		--readFilesCommand ${params.config.alignment.readFilesCommand} \
		--outReadsUnmapped ${params.config.alignment.outReadsUnmapped} \
		--outSAMunmapped  ${params.config.alignment.outSAMunmapped} \
		--outSAMtype ${params.config.alignment.outSAMtype} \
		--outSAMattrRGline ID:$sample_name.$lane_id LB:library PL:Illumina SM:$sample_name \
		--limitBAMsortRAM ${params.config.memory.alignmentBAMSort} \
		--runThreadN ${params.config.alignment.runThreadN} \
		--chimSegmentMin ${params.config.alignment.chimSegmentMin}  \
		--chimJunctionOverhangMin ${params.config.alignment.chimJunctionOverhangMin} \
		--chimScoreMin ${params.config.alignment.chimScoreMin} \
		--chimScoreDropMax ${params.config.alignment.chimScoreDropMax} \
		--chimScoreSeparation ${params.config.alignment.chimScoreSeparation} \
		--chimScoreJunctionNonGTAG ${params.config.alignment.chimScoreJunctionNonGTAG} \
		--sjdbOverhang ${params.config.alignment.sjdbOverhang} \
		--seedSearchStartLmax ${params.config.alignment.seedSearchStartLmax} \
		--outFileNamePrefix ${lane_id}. \
		--outTmpDir ${params.javatmp}/${params.uuid}/STAR-${runid}-${sample_name}.${lane_id}
		
	ev=\$?

	sleep 60

} 2>&1

script_exit
