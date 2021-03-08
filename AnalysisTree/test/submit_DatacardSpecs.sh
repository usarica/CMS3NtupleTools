#!/bin/bash

date=$1
period=$2
templateVersion=$3

script=produceDatacardSpecs.cc
function=runDatacardChain
jobdate="${date}_DCSpecs"
arguments='"<period>","<templateVersion>","<strdate>",<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<templateVersion>/$templateVersion}"

for includeBoostedHadVHCategory in false true; do
for includeResolvedHadVHCategory in false true; do
  strargs="${arguments}"
  strargs="${strargs/<includeBoostedHadVHCategory>/${includeBoostedHadVHCategory}}"
  strargs="${strargs/<includeResolvedHadVHCategory>/${includeResolvedHadVHCategory}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
done
done
