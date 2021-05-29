#!/bin/bash

date=$1
period=$2
ntupleVersion=$3

script=produceSimBkgTemplates.cc
function=runTemplateChain
jobdate="${date}_SimBkgTpls"
arguments='"<period>","<ntupleVersion>","<strdate>",<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>,"<strSampleSet_target>"'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"

for includeBoostedHadVHCategory in false true; do
for includeResolvedHadVHCategory in false; do
for strSampleSet_target in qqZZ_offshell qqWZ_offshell tZX; do
  strargs="${arguments}"
  strargs="${strargs/<includeBoostedHadVHCategory>/${includeBoostedHadVHCategory}}"
  strargs="${strargs/<includeResolvedHadVHCategory>/${includeResolvedHadVHCategory}}"
  strargs="${strargs/<strSampleSet_target>/${strSampleSet_target}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
done
done
done
