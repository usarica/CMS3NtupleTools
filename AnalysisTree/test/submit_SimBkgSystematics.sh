#!/bin/bash

date=$1
period=$2
prodVersion=$3

script=produceSimBkgEstimates.cc
function=produceSystematicsReweighting
jobdate="${date}_SimBkg_SystematicsReweighting"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>",<theGlobalSyst>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<prodVersion>/$prodVersion}"

declare -a MCSysts=( tPythiaScaleDn tPythiaScaleUp )
declare -a SimSamples=()
if [[ "$period" == "2017" ]]; then
  SimSamples=( "qqZZ_2l2nu_mZ_18-inf" "qqZZ_2l2nu" "qqZZ_4l" "qqZZ_4l_ext" "qqWZ_3lnu_POWHEG_mll_0p1-inf" "qqWZ_3lnu_POWHEG" "qqWW_lnu2q" "qqWW_lnu2q_ext" "qqWW_2l2nu" "qqWW_2l2nu_ext" )
elif [[ "$period" == "2018" ]]; then
  SimSamples=( "qqZZ_2l2nu_mZ_18-inf" "qqZZ_2l2nu" "qqZZ_2l2nu_ext" "qqZZ_4l" "qqWZ_3lnu_POWHEG_mll_0p1-inf" "qqWZ_3lnu_POWHEG" "qqWW_lnu2q" "qqWW_2l2nu" )
fi
for sample in "${SimSamples[@]}"; do
  for syst in "${MCSysts[@]}"; do
    strargs="${arguments}"
    strargs="${strargs/<strSampleSet>/${sample}}"
    strargs="${strargs/<theGlobalSyst>/${syst}}"

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
  done
done
