#!/bin/bash

date=$1
period=$2
ntupleVersion=$3

script=produceNRBTemplates.cc
function=runTemplateChain
jobdate="${date}_NRBTpls"
arguments='"<period>","<ntupleVersion>","<strdate>",<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"

for includeBoostedHadVHCategory in false true; do
for includeResolvedHadVHCategory in false true; do
  strargs="${arguments}"
  strargs="${strargs/<includeBoostedHadVHCategory>/${includeBoostedHadVHCategory}}"
  strargs="${strargs/<includeResolvedHadVHCategory>/${includeResolvedHadVHCategory}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
done
done
