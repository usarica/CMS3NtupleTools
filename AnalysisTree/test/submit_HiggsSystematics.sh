#!/bin/bash

date=$1
period=$2
prodVersion=$3

script=produceHiggsEstimates.cc
function=produceSystematicsReweighting
jobdate="${date}_SystematicsReweighting"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>",<theGlobalSyst>'
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

declare -a MCSysts=( )
declare -a SimSamples=( $(echo $strSampleGroup) )
for sample in "${SimSamples[@]}"; do
  MCSysts=( tPythiaTuneDn tPythiaTuneUp )
  if [[ "$period" == "2016" ]]; then
    MCSysts+=( tPythiaScaleDn tPythiaScaleUp )
  else
    if [[ "$sample" == "GGH"* ]] || [[ "$sample" == "VBF"* ]]; then
      MCSysts+=( tAsMZDn tAsMZUp )
    fi
  fi
  if [[ "$sample" == "GGH"* ]]; then
    MCSysts+=( tHardJetsDn tHardJetsUp tPDFScaleDn tPDFScaleUp tQCDScaleDn tQCDScaleUp tPDFReplicaDn tPDFReplicaUp )
  fi

  for syst in "${MCSysts[@]}"; do
    strargs="${arguments}"
    strargs="${strargs/<strSampleSet>/${sample}}"
    strargs="${strargs/<theGlobalSyst>/${syst}}"

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
  done
done
