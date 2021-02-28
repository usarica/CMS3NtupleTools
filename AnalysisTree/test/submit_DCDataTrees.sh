#!/bin/bash

date=$1
period=$2
prodVersion=$3
ntupleVersion=$4

script=produceDCDataTrees.cc
function=runDistributionsChain
jobdate="${date}_DCDataTrees"
arguments='"<period>","<prodVersion>","<ntupleVersion>","<strdate>",<AChypo>,<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<prodVersion>/$prodVersion}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"

for includeBoostedHadVHCategory in false true; do
for includeResolvedHadVHCategory in false true; do
for hypo in kSM kA3 kA2 kL1 kL1ZGs; do
  strargs="${arguments}"
  strargs="${strargs/<AChypo>/${hypo}}"
  strargs="${strargs/<includeBoostedHadVHCategory>/${includeBoostedHadVHCategory}}"
  strargs="${strargs/<includeResolvedHadVHCategory>/${includeResolvedHadVHCategory}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
done
done
done