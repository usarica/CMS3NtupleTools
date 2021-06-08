#!/bin/bash

date=$1
period=$2
ntupleVersion=$3
channel=""
if [[ "$4" != "" ]]; then
  channel="$4"
else
  echo "A channel is mandatory. Please pass ZZ or WZ."
  exit 1
fi
if [[ "$channel" == "WZ"* ]] || [[ "$channel" == "ZW"* ]]; then
  channel="ZWTo3L1Nu"
else
  channel="ZZTo2L2Nu"
fi

script=produceHiggsTemplates.cc
function=runTemplateChain_${channel}
jobdate="${date}_HiggsTpls_${channel}"
arguments='"<strSampleSet>","<period>","<ntupleVersion>","<strdate>",<AChypo>,<syst>,<includeBoostedHadVHCategory>,<includeResolvedHadVHCategory>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"

declare -a procs=( ggZZ_offshell VVVV_offshell )
if [[ "$channel" == "ZW"* ]]; then
  procs=( VVVV_offshell )
fi

declare -a systs=()

for proc in "${procs[@]}"; do
systs=(
    sNominal \
    tPDFScaleDn tPDFScaleUp \
    tQCDScaleDn tQCDScaleUp \
    tAsMZDn tAsMZUp \
    tPDFReplicaDn tPDFReplicaUp \
    tPythiaScaleDn tPythiaScaleUp \
    tPythiaTuneDn tPythiaTuneUp \
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

if [[ "$channel" == "ZW"* ]] && [[ "$includeBoostedHadVHCategory" == "true" ]]; then
  continue
fi

for includeResolvedHadVHCategory in false; do

if [[ "$channel" == "ZW"* ]] && [[ "$includeResolvedHadVHCategory" == "true" ]]; then
  continue
fi

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
