#!/bin/bash

date=$1
period=$2
prodVersion=$3
script=produceGJetsMETResolution.cc
function=getTrees
jobdate="${date}_METResolution"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>",<ichunk>,<nchunks>,<theGlobalSyst>'
arguments="${arguments/<strdate>/${date}}"
arguments="${arguments/<prodVersion>/$prodVersion}"

declare -a dataPeriods=( )
if [[ "$period" == "2018" ]]; then
  dataPeriods=( 2018A 2018B 2018C 2018D )
elif [[ "$period" == "2017" ]]; then
  dataPeriods=( 2017B 2017C 2017D 2017E 2017F )
elif [[ "$period" == "2016" ]]; then
  dataPeriods=( 2016B 2016C 2016D 2016E 2016F 2016G 2016H )
fi

csvfile="skimSamples_${period}.csv"
for sample in $(readCMS3SkimSamplesFromCSV.py --csv=${csvfile} --sim --tree_req="SinglePhoton"); do
  if [[ "$sample" == *"WJetsToLNu_HT"* ]]; then
    continue
  fi

  sampleDir=${sample//MINIAODSIM}

  echo "====="
  skimdir="/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/${prodVersion}/${sampleDir}"
  if [[ ! -d $skimdir ]]; then
    echo "$skimdir does not exist"
    continue
  fi

  let nfiles=$(ls $skimdir | grep -e ".root" | wc -l)
  echo "$sampleDir has $nfiles files"

  nchunks=$((nfiles / 4))
  nfeff=$((nchunks * 4))
  if [[ $nfeff -lt $nfiles ]];then
    nchunks=$((nchunks + 1))
  fi
  echo " - There will be $nchunks chunks"
  if [[ $nchunks -eq 0 ]]; then
    echo "Oh-uh..."
  fi

  for dataperiod in "${dataPeriods[@]}"; do
  for syst in sNominal eJECDn eJECUp eJERDn eJERUp ePUDn ePUUp; do

  let ichunk=0
  while [[ $ichunk -lt $nchunks ]]; do
    strargs="${arguments}"
    strargs="${strargs/<strSampleSet>/${sample}}"
    strargs="${strargs/<ichunk>/${ichunk}}"
    strargs="${strargs/<nchunks>/${nchunks}}"
    strargs="${strargs/<theGlobalSyst>/${syst}}"
    strargs="${strargs/<period>/$dataperiod}"

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" job_flavor="workday" no-proxycheck
    let ichunk=$ichunk+1
  done

  done
  done
  
  echo "====="
done

for dataperiod in "${dataPeriods[@]}"; do
  sample="Run${dataperiod}"

  strargs="${arguments}"
  strargs="${strargs/<strSampleSet>/${sample}}"
  strargs="${strargs/<ichunk>/0}"
  strargs="${strargs/<nchunks>/0}"
  strargs="${strargs/<theGlobalSyst>/sNominal}"
  strargs="${strargs/<period>/$dataperiod}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" job_flavor="workday" no-proxycheck
done
