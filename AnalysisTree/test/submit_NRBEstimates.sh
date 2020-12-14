#!/bin/bash

date=$1
period=$2
prodVersion=$3
script=produceNRBEstimates.cc
function=runDistributionsChain
jobdate="NRBEstimates_${date}"
arguments='"<period>","<prodVersion>","<strdate>",<theGlobalSyst>'
arguments="${arguments/<period>/$period}"
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<prodVersion>/$prodVersion}"

# Run only on those that change data and MC alike
declare -a systs=( \
  sNominal \
  eEleEffStatDn eEleEffStatUp \
  eEleEffSystDn eEleEffSystUp \
  eEleEffAltMCDn eEleEffAltMCUp \
  eMuEffStatDn eMuEffStatUp \
  eMuEffSystDn eMuEffSystUp \
  eMuEffAltMCDn eMuEffAltMCUp \
  eTriggerEffDn eTriggerEffUp \
)


for syst in "${systs[@]}"; do
  strargs="${arguments}"
  strargs="${strargs/<theGlobalSyst>/${syst}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
done
