#!/bin/sh

exec < /uscms_data/d1/jmuelmen/Zjets-madgraph.txt
exec > /uscms_data/d1/jmuelmen/Zjets-madgraph-files.txt

declare -a line
while read -a line ; do
    lumi=${line[2]} 
    file=`$DBSCMD_HOME/dbsCommandLine.py -c search --query="find file where lumi = $lumi and dataset like /ZJets-madgraph/Summer08_IDEAL_V11_redigi_v1/GEN-SIM-RECO" | grep /store`
    echo ${line[0]} ${line[1]} ${line[2]} $file
done
