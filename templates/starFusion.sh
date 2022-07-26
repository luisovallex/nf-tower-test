#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

echo "chimeric sam file: ${chimeric_sam}, jxns: ${chimeric_junctions_merged}"
echo "genome dir: ${fusiongenome}"

{

	$params.starfusion_path/STAR-Fusion \
		--genome_lib_dir ${fusiongenome} \
		-J ${chimeric_junctions_merged} \
		${params.config.starFusion.starFusionArgs} \
		-O fusion_results
	ev=\$?
	sleep 60

	

} 2>&1

script_exit
