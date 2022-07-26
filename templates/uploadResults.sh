#!/bin/bash

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{
	#echo "XXX diabled SQS"

	#echo "uploading final results from ${runid}/${sample_name} to ${params.output}..."
	#aws s3 sync ${params.results}/${runid}/${sample_name}/ ${params.output}/${runid}/${sample_name}
	echo "output files can be found in ${params.results}/${runid}/${sample_name}"

	#echo "sending completion message to SQS..."
	#aws sqs send-message  --queue-url ${params.sqs_queue} --message-body "PipelineComplete:${runid}:${sample_name}"

	echo "syncing logs..."
	#aws s3 sync ${params.logs} ${params.output}/logs/`date "+%Y%m%d"`-`hostname --ip-address`-`hostname`/
	aws s3 sync ${params.logs} ${params.output}/logs/`date "+%Y%m%d"`-`hostname`/

	#aws sqs send-message  --queue-url ${params.sqs_queue} --message-body "LogSyncComplete:${runid}:${sample_name}"-`hostname`
	echo "complete!"
	ev=\$?

	

} 2>&1

script_exit
