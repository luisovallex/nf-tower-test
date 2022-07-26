#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

echo "${ds_percent}"

{

	/usr/local/pipeline/packages/samtools-1.0/samtools view -b -s ${ds_percent} ${bam} > downsampled.8m.bam

	java ${params.config.memory.markDuplicatesStart} ${params.config.memory.markDuplicatesMax} -Djava.io.tmpdir=$params.javatmp -jar $params.picard_path/MarkDuplicates.jar \
			INPUT=downsampled.8m.bam \
			OUTPUT=dedupped.8m.bam \
			METRICS_FILE=${params.config.downsampling8m.METRICS_FILE} \
			CREATE_INDEX=${params.config.downsampling8m.CREATE_INDEX} \
			VALIDATION_STRINGENCY=${params.config.downsampling8m.VALIDATION_STRINGENCY}
	ev=\$?
	sleep 60
	

} 2>&1

script_exit
