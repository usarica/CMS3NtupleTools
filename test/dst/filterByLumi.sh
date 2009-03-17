#!/bin/sh

exec < Zjets-madgraph-files.txt

declare -a line
while read -a line ; do
    export RUN=${line[0]}
    export EVENT=${line[1]}
    export LUMI=${line[2]} 
    nfiles=$((${#line[*]} - 3))
    while [ $nfiles -gt 0 ] ; do
	export INFILE=${line[$((nfiles + 2))]}
	export OUTFILE=file:LUMI-${LUMI}-`basename $INFILE`
	echo stripping run/lumi section $RUN/$LUMI from $INFILE to $OUTFILE
	( cmsRun filterByLumi.py ) &
	nfiles=$(($nfiles - 1))
    done 
done

wait
echo all done 
