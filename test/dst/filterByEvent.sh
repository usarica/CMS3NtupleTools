#!/bin/sh

exec < Zjets-madgraph-files.txt

declare -a line
declare -a outputs
while read -a line ; do
    export RUN=${line[0]}
    export EVENT=${line[1]}
    export LUMI=${line[2]} 
    nfiles=$((${#line[*]} - 3))
    while [ $nfiles -gt 0 ] ; do
	export INFILE=${line[$((nfiles + 2))]}
	export OUTFILE=file:LUMI-${LUMI}-`basename $INFILE`
	export OUTFILE2=${OUTFILE/LUMI/RUN-${RUN}-EVENT-${EVENT}-LUMI}
	echo stripping run/event $RUN/$EVENT from $OUTFILE to $OUTFILE2
	( INFILE=$OUTFILE 
	    OUTFILE=$OUTFILE2
	    if [ -e ${INFILE/file:/} ] ; then
		cmsRun filterByEvent.py 
		rm ${INFILE/file:/} 
	    fi
	) &
	outputs[${#outputs[*]}]=${OUTFILE2}
	nfiles=$(($nfiles - 1))
    done 
done

wait
echo all done, outputs: 
iout=0
while [ $iout -lt ${#outputs[*]} ] ; do
    if [ -e ${outputs[$iout]/file:/} ] ; then
	echo -n \"${outputs[$iout]}\",
    fi
    iout=$(($iout + 1))
done
echo
