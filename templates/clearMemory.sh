#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} $runid

{

	ipcs -m
	ev=\$?

	

} 2>&1


script_exit
