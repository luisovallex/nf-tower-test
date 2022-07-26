#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$runid"

{

	echo "startup for runid $runid"
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
