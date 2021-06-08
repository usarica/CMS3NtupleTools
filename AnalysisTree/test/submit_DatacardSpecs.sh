#!/bin/bash

date=$1
period=$2
templateVersion=$3
channel=""
if [[ "$4" == "ZZ"* ]]; then
  channel="ZZTo2L2Nu"
elif [[ "$4" == "ZW"* ]] || [[ "$4" == "WZ"* ]]; then
  channel="ZWTo3L1Nu"
fi
if [[ "${channel}" == "" ]]; then
  echo "Channel string is empty."
  exit 1
fi

script=produceDatacardSpecs_${channel}.cc
function=runDatacardChain
jobdate="${date}_DCSpecs_${channel}"
arguments='"<period>","<templateVersion>","<strdate>",<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<templateVersion>/$templateVersion}"

for includeBoostedHadVHCategory in false true; do
  if [[ "${channel}" == "ZW"* ]] && [[ "${includeBoostedHadVHCategory}" == "true" ]]; then
    continue
  fi

  for includeResolvedHadVHCategory in false; do
    if [[ "${channel}" == "ZW"* ]] && [[ "${includeResolvedHadVHCategory}" == "true" ]]; then
      continue
    fi

    strargs="${arguments}"
    strargs="${strargs/<includeBoostedHadVHCategory>/${includeBoostedHadVHCategory}}"
    strargs="${strargs/<includeResolvedHadVHCategory>/${includeResolvedHadVHCategory}}"

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
  done
done
