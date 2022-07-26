#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	$params.samtools_path/samtools sort \
		${params.config.sortByName.samtoolsSortArgs} ${params.config.memory.samtoolsSort} \
		${bam} \
		${params.config.sortByName.samtoolsSortType}
	ev=\$?

	sleep 90

} 2>&1

script_exit
