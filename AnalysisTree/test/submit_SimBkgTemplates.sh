#!/bin/bash

date=$1
period=$2
ntupleVersion=$3
strSampleSet_target=""
for arg in "$@"; do
  if [[ "$arg" == "do_only="* ]]; then
    strSampleSet_target=${arg#*=}
  fi
done


script=produceSimBkgTemplates.cc
function=runTemplateChain
jobdate="${date}_SimBkgTpls"
arguments='"<period>","<ntupleVersion>","<strdate>",<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>,"<strSampleSet_target>"'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"
arguments="${arguments/<strSampleSet_target>/$strSampleSet_target}"

for includeBoostedHadVHCategory in false true; do
for includeResolvedHadVHCategory in false true; do
  strargs="${arguments}"
  strargs="${strargs/<includeBoostedHadVHCategory>/${includeBoostedHadVHCategory}}"
  strargs="${strargs/<includeResolvedHadVHCategory>/${includeResolvedHadVHCategory}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
done
done
