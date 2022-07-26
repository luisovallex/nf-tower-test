#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	java ${params.config.memory.markDuplicatesStart} ${params.config.memory.markDuplicatesMax} -Djava.io.tmpdir=$params.javatmp -jar $params.picard_path/MarkDuplicates.jar \
		INPUT=${bam} \
		OUTPUT=dedupped.bam \
		METRICS_FILE=${params.config.mark_duplicates.METRICS_FILE} \
		CREATE_INDEX=${params.config.mark_duplicates.CREATE_INDEX} \
		VALIDATION_STRINGENCY=${params.config.mark_duplicates.VALIDATION_STRINGENCY}
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
