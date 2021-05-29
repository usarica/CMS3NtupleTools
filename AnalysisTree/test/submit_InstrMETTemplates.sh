#!/bin/bash

date=$1
period=$2
ntupleVersion=$3
skipIntermediates="false"
for arg in "$@"; do
  if [[ "$arg" == "skipIntermediates="* ]]; then
    skipIntermediates=${arg#*=}
  fi
done
if [[ "${skipIntermediates}" != "true" ]] && [[ "${skipIntermediates}" != "false" ]]; then
  echo "skipIntermediates must be 'true' or 'false'"
fi


script=produceInstrMETTemplates.cc
function=runTemplateChain
jobdate="${date}_InstrMETTpls"
arguments='"<period>","<ntupleVersion>","<strdate>",<AChypo>,<dilepton_id>,<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>,<skipIntermediates>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"
arguments="${arguments/<skipIntermediates>/$skipIntermediates}"

for includeBoostedHadVHCategory in false true; do
for includeResolvedHadVHCategory in false; do
for AChypo in kSM kA3 kA2 kL1 kL1ZGs; do
for dilepton_id in -121 -169; do
  strargs="${arguments}"
  strargs="${strargs/<AChypo>/${AChypo}}"
  strargs="${strargs/<dilepton_id>/${dilepton_id}}"
  strargs="${strargs/<includeBoostedHadVHCategory>/${includeBoostedHadVHCategory}}"
  strargs="${strargs/<includeResolvedHadVHCategory>/${includeResolvedHadVHCategory}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
done
done
done
done
