#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{

	$params.samtools_path/samtools view \
		-bh -F 0x10 ${bam} > nomult.bam

	ls -alt nomult.bam

	$params.samtools_path/samtools index \
		nomult.bam

	# rename index file
	mv nomult.bam.bai nomult.bai

	# XXX not parallel, for now

	# index stats
	$params.samtools_path/samtools idxstats \
		${bam} > dedupped.stats.txt

	$params.samtools_path/samtools idxstats \
		nomult.bam > nomult.stats.txt


	# stats
	$params.samtools_path/samtools stats \
		${bam} > dedupped.stattbl.txt

	$params.samtools_path/samtools stats \
		nomult.bam > nomult.stattbl.txt


	# stats, remove duplicates
	$params.samtools_path/samtools stats --remove-dups \
		${bam} > dedupped.stattbl2.txt

	$params.samtools_path/samtools stats --remove-dups \
		nomult.bam > nomult.stattbl2.txt


	# flagstat
	$params.samtools_path/samtools flagstat \
		${bam} > dedupped.flags.txt

	$params.samtools_path/samtools flagstat \
		nomult.bam > nomult.flags.txt


	# header
	$params.samtools_path/samtools view -h \
		nomult.bam 1:10000000-11000000 > nomult.header.txt
	ev=\$?

	sleep 60

	

} 2>&1

script_exit
