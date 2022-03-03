#!/bin/bash

date=$1
period=$2
prodVersion=$3
readFromFile=""
if [[ $# -gt 3 ]]; then
  readFromFile=$4
  echo "Reading samples from ${readFromFile=}"
fi
script=scanTrees.cc
function=scanTrees
jobdate="ValidTreeScan_${date}"
arguments='"<strSampleSet>","<period>","<prodVersion>",true,true'
arguments="${arguments/<period>/$period}"
arguments="${arguments/<prodVersion>/$prodVersion}"

proddir="/ceph/cms/store/user/usarica/Offshell_2L2Nu/Production/${prodVersion}"
proddir=$( ExecuteCompiledCommand GetStandardHostPathToStore ${proddir} t2.ucsd.edu )
if [[ "$( ExecuteCompiledCommand DirectoryExists ${proddir} )" == "false" ]]; then
  echo "$proddir does not exist"
  exit 1
fi

declare -a checklist=( )

if [[ "${readFromFile}" != "" ]]; then
  while IFS='' read -r line || [[ -n "$line" ]]; do
    if [[ "${line}" == "" ]]; then
      continue
    fi
    checklist+=( ${line} )
  done < "${readFromFile}"
else
  for sdir in $( ExecuteCompiledCommand lsdir ${proddir} ); do
    for ssdir in $( ExecuteCompiledCommand lsdir ${proddir}/${sdir} ); do
      strSampleSet="/${sdir}/${ssdir}"
      if [[ "${ssdir}" == *"RunII"* ]]; then
        strSampleSet="${strSampleSet}/MINIAODSIM"
      else
        strSampleSet="${strSampleSet}/MINIAOD"
      fi
      checklist+=( ${strSampleSet} )
    done
  done
fi



for strSampleSet in "${checklist[@]}"; do
  echo "Found sample set ${strSampleSet}"

  strargs="${arguments}"
  strargs="${strargs/<strSampleSet>/${strSampleSet}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
done
