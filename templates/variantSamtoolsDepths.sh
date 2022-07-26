#!/bin/bash

# depths: operates on calibrated BAM, generates nomult.bam

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	$params.samtools_path/samtools view -bh -F 0x100 \
		${bam} > nomult.bam

	$params.samtools_path/samtools index nomult.bam

	$params.samtools_path/samtools depth \
		-b ${resources}/bambooTCGA.ThyroseqV2V2Rasal854sites30Jul2015.bed \
		nomult.bam > ThyroVIPv2.nomult.depths.txt

	$params.samtools_path/samtools depth \
		-b ${resources}/bambooTCGA.ThyroseqV2V2Rasal854sites30Jul2015.bed \
		${bam} > ThyroVIPv2.depths.txt
	ev=\$?
	sleep 60
	

} 2>&1

script_exit
