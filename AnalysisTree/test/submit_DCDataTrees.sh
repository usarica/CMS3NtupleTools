#!/bin/bash

date=$1
period=$2
prodVersion=$3
ntupleVersion=$4
channel=$5

arguments=""
if [[ "$channel" == "WZ"* ]] || [[ "$channel" == "ZW"* ]]; then
  channel="ZWTo3L1Nu"
  arguments='"<period>","<prodVersion>","<ntupleVersion>","<strdate>",<lepid>,<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>'
else
  channel="ZZTo2L2Nu"
  arguments='"<period>","<prodVersion>","<ntupleVersion>","<strdate>",<AChypo>,<lepid>,<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>'
fi

script=produceDCDataTrees.cc
function=getTrees_${channel}
jobdate="${date}_DCDataTrees_${channel}"

arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<prodVersion>/$prodVersion}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"

lepids=()
if [[ "$channel" == "ZZTo2L2Nu" ]]; then
  lepids=( -121 -169 )
else
  lepids=( -11 -13 )
fi

for lepid in "${lepids[@]}"; do
for includeBoostedHadVHCategory in false true; do
for includeResolvedHadVHCategory in false; do
for hypo in kSM kA3 kA2 kL1 kL1ZGs; do
  # There is no need for an AC hypothesis in ZW data since discriminants are not used.
  if [[ "$channel" == "ZWTo3L1Nu" ]] && [[ "$hypo" != "kSM" ]]; then
    continue
  fi

  strargs="${arguments}"
  strargs="${strargs/<AChypo>/${hypo}}"
  strargs="${strargs/<lepid>/${lepid}}"
  strargs="${strargs/<includeBoostedHadVHCategory>/${includeBoostedHadVHCategory}}"
  strargs="${strargs/<includeResolvedHadVHCategory>/${includeResolvedHadVHCategory}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
done
done
done
done
