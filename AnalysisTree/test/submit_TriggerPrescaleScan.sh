#!/bin/bash

date=$1
prodVersion=$2
script=scanTriggerPrescales.cc
function=scanTriggerPrescales
jobdate="${date}_TriggerPrescaleScan"
arguments='"<period>","<prodVersion>"'


for period in 2016 2017 2018; do
  prodtag=${prodVersion}_${period}

  declare -a dataPeriods

  if [[ "$period" == "2018" ]]; then
    dataPeriods=( 2018A 2018B 2018C 2018D )
  elif [[ "$period" == "2017" ]]; then
    dataPeriods=( 2017B 2017C 2017D 2017E 2017F )
  elif [[ "$period" == "2016" ]]; then
    dataPeriods=( 2016B 2016C 2016D 2016E 2016F 2016G 2016H )
  fi


  for dataperiod in "${dataPeriods[@]}"; do
    strargs="${arguments}"
    strargs="${strargs/<prodVersion>/$prodtag}"
    strargs="${strargs/<period>/$dataperiod}"

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
  done
done
