#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	java ${params.config.memory.mergeStart} ${params.config.memory.mergeMax} -Djava.io.tmpdir=$params.javatmp -jar $params.picard_path/MergeSamFiles.jar \
		${bam.collect({'INPUT=' + it + ' '}).join()} \
		OUTPUT=merged.bam \
		USE_THREADING=${params.config.merge.USE_THREADING} \
                CREATE_INDEX=${params.config.merge.CREATE_INDEX} \
		VALIDATION_STRINGENCY=${params.config.merge.VALIDATION_STRINGENCY}

	## generate index bam
	#/usr/local/pipeline/packages/samtools-1.0/samtools index merged.bam

	ev=\$?

	touch cleanup
	sleep 60

	

} 2>&1

script_exit
