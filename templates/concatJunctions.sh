#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	echo "concatJunctions step, chimeric_junctions"
	/usr/local/setup/nextflow/bin/concat_junctions.pl ${chimeric_junctions.collect({it + ' '}).join()} >> ${sample_name}.chimeric_junctions_merged_stdout.txt
	ev=\$?

	ls -alt chimeric_junctions_merged

	sleep 90

	

} 2>&1

script_exit
