#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

echo "${task.process}: input was ${chimeric_sam}" 

{

	java ${params.config.memory.mergeChimericStart} ${params.config.memory.mergeChimericMax} -Djava.io.tmpdir=$params.javatmp -jar $params.picard_path/MergeSamFiles.jar \
		${chimeric_sam.collect({'INPUT=' + it + ' '}).join()} \
		OUTPUT=concat.chimeric.out.sam \
		USE_THREADING=${params.config.mergeChimericSams.USE_THREADING} \
		VALIDATION_STRINGENCY=${params.config.mergeChimericSams.VALIDATION_STRINGENCY}
	ev=\$?

	sleep 90

} 2>&1

script_exit
