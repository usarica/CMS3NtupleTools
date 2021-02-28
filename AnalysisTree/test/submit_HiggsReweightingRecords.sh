#!/bin/bash

date=$1
period=$2
prodVersion=$3

script=produceHiggsEstimates.cc
function=produceReweightingRecords
jobdate="${date}_ReweightingRecords"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>"'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<prodVersion>/$prodVersion}"

strSampleGroup=""
if [[ "$period" == "2018"* ]]; then
  strSampleGroup="${strSampleGroup} GGH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} GGH_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2Nu2X_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2L2Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo4Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_WWTo2L2Nu_POWHEG"
  #strSampleGroup="${strSampleGroup} ZH_HToLNuQQ_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Nu_POWHEG"
  #strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_HToWW_2LOSFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Nu_POWHEG"
  #strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_HToWW_2LOSFilter_POWHEG"
elif [[ "$period" == "2017"* ]]; then
  strSampleGroup="${strSampleGroup} GGH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} GGH_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2Nu2X_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2L2Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo4Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_WWTo2L2Nu_POWHEG"
  #strSampleGroup="${strSampleGroup} ZH_HToLNuQQ_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Nu_POWHEG"
  #strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_HToWW_2LOSFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Nu_POWHEG"
  #strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_HToWW_2LOSFilter_POWHEG"
elif [[ "$period" == "2016"* ]]; then
  strSampleGroup="${strSampleGroup} GGH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} GGH_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2Nu2X_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2L2Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo4Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_WWTo2L2Nu_POWHEG"
  #strSampleGroup="${strSampleGroup} ZH_HToLNuQQ_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Nu_POWHEG"
  #strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_HToWW_2LOSFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Nu_POWHEG"
  #strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_HToWW_2LOSFilter_POWHEG"
fi

declare -i SimSamples=( $(echo $strSampleGroup) )
for sample in "${SimSamples[@]}"; do
  strargs="${arguments}"
  strargs="${strargs/<strSampleSet>/${sample}}"
  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
done
