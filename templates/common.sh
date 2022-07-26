# bash
# library of common function for Nextflow template scripts

test_echo () {
	date
	echo $1
}

begin_timing () {
	begin=$(date +%s)
}

end_timing () {
	end=$(date +%s)
	timediff=$(( $end - $begin ))
	printf "elapsed time:\t%d seconds\n" $timediff
}	

script_begin() {
	# passed task name and sample name
	echo ""
	echo "----"
	printf "%s: %s\n" $1 $2
	echo "runhost is: `hostname -s`"
	begin_timing
}

script_exit() {
	end_timing
	sleep 5

	if ! [[ "$ev" =~ ^[0-9]+$ ]]
		then
			echo "exit value variable improperly set, or missing!"
			exit -1
	fi

	printf "exit value is:\t%d\n" $ev
	echo ""
	exit $ev
}

# cat /dev/notthere
# ev=\$?
# sleep 1
# exit \$ev


