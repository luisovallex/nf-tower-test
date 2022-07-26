#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	echo "removing result files...."
	#/bin/rm -rf ${params.results}/${sample_name}/
	echo "would remove ${params.results}/${runid}/${sample_name}/"
	echo "complete!"
	ev=\$?

	sleep 60

} 2>&1

script_exit
