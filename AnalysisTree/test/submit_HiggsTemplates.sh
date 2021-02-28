#!/bin/bash

date=$1
period=$2
ntupleVersion=$3

script=produceHiggsTemplates.cc
function=runTemplateChain
jobdate="${date}_HiggsTpls"
arguments='"<strSampleSet>","<period>","<ntupleVersion>","<strdate>",<AChypo>,<syst>,<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"

declare -a procs=( ggZZ_offshell VVVV_offshell )
declare -a systs=()

for proc in "${procs[@]}"; do
systs=(
    sNominal \
    tPDFScaleDn tPDFScaleUp \
    tQCDScaleDn tQCDScaleUp \
    tAsMZDn tAsMZUp \
    tPDFReplicaDn tPDFReplicaUp \
    tPythiaScaleDn tPythiaScaleUp \
    ePhoEffDn ePhoEffUp \
    eMETDn eMETUp \
    eJECDn eJECUp \
    eJERDn eJERUp \
    ePUDn ePUUp \
    ePUJetIdEffDn ePUJetIdEffUp \
    eBTagSFDn eBTagSFUp \
    eTriggerEffDn eTriggerEffUp \
    eEleEffStatDn eEleEffStatUp \
    eEleEffSystDn eEleEffSystUp \
    eEleEffAltMCDn eEleEffAltMCUp \
    eMuEffStatDn eMuEffStatUp \
    eMuEffSystDn eMuEffSystUp \
    eMuEffAltMCDn eMuEffAltMCUp \
)
if [[ "${period}" == "2016"* ]] || [[ "${period}" == "2017"* ]]; then
  systs+=( eL1PrefiringDn eL1PrefiringUp )
fi
if [[ "${proc}" == "gg"* ]]; then
  systs+=( tHardJetsDn tHardJetsUp )
fi

for syst in "${systs[@]}"; do
for includeBoostedHadVHCategory in false true; do
for includeResolvedHadVHCategory in false true; do
for hypo in kSM kA3 kA2 kL1 kL1ZGs; do
  strargs="${arguments}"
  strargs="${strargs/<strSampleSet>/${proc}}"
  strargs="${strargs/<syst>/${syst}}"
  strargs="${strargs/<AChypo>/${hypo}}"
  strargs="${strargs/<includeBoostedHadVHCategory>/${includeBoostedHadVHCategory}}"
  strargs="${strargs/<includeResolvedHadVHCategory>/${includeResolvedHadVHCategory}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
done
done
done
done

done