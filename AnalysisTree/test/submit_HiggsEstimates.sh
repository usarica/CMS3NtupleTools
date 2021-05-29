#!/bin/bash

date=$1
period=$2
prodVersion=$3
ntupleVersion=$4
rewgtRcdVersion=$5

script=produceHiggsEstimates.cc
function=runDistributionsChain
jobdate="${date}_HiggsEstimates"
arguments='"<strSampleSet>","<period>","<prodVersion>","<ntupleVersion>","<rewgtRcdVersion>","<strdate>",<theGlobalSyst>,<applyPUIdToAK4Jets>,<applyTightLeptonVetoIdToAK4Jets>,<use_MET_Puppi>,<use_MET_XYCorr>,<use_MET_JERCorr>,<use_MET_ParticleMomCorr>,<use_MET_p4Preservation>,<use_MET_corrections>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<period>/$period}"
arguments="${arguments/<prodVersion>/$prodVersion}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"
arguments="${arguments/<rewgtRcdVersion>/$rewgtRcdVersion}"
arguments="${arguments/<applyPUIdToAK4Jets>/true}"
arguments="${arguments/<applyTightLeptonVetoIdToAK4Jets>/false}"
arguments="${arguments/<use_MET_Puppi>/false}"
arguments="${arguments/<use_MET_XYCorr>/true}"
arguments="${arguments/<use_MET_JERCorr>/false}"
arguments="${arguments/<use_MET_ParticleMomCorr>/true}"
arguments="${arguments/<use_MET_p4Preservation>/true}"
arguments="${arguments/<use_MET_corrections>/true}"

declare -a systs=( \
    sNominal \
    tPDFScaleDn tPDFScaleUp \
    tQCDScaleDn tQCDScaleUp \
    tAsMZDn tAsMZUp \
    tPDFReplicaDn tPDFReplicaUp \
    tPythiaScaleDn tPythiaScaleUp \
    tPythiaTuneDn tPythiaTuneUp \
    tHardJetsDn tHardJetsUp \
    eEleEffStatDn eEleEffStatUp \
    eEleEffSystDn eEleEffSystUp \
    eEleEffAltMCDn eEleEffAltMCUp \
    eMuEffStatDn eMuEffStatUp \
    eMuEffSystDn eMuEffSystUp \
    eMuEffAltMCDn eMuEffAltMCUp \
    ePhoEffDn ePhoEffUp \
    eMETDn eMETUp \
    eJECDn eJECUp \
    eJERDn eJERUp \
    ePUDn ePUUp \
    ePUJetIdEffDn ePUJetIdEffUp \
    eBTagSFDn eBTagSFUp \
    eTriggerEffDn eTriggerEffUp \
)
if [[ "$period" == "2016"* ]] || [[ "$period" == "2017"* ]]; then
  systs+=( eL1PrefiringDn eL1PrefiringUp )
fi


strSampleGroup=""
if [[ "$period" == "2018"* ]]; then
  strSampleGroup="${strSampleGroup} GGH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} GGH_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2Nu2X_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2L2Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo4Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HToLNuQQ_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_HToWW_2LOSFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_HToWW_2LOSFilter_POWHEG"
elif [[ "$period" == "2017"* ]]; then
  strSampleGroup="${strSampleGroup} GGH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} GGH_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2Nu2X_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2L2Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo4Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HToLNuQQ_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_HToWW_2LOSFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_HToWW_2LOSFilter_POWHEG"
elif [[ "$period" == "2016"* ]]; then
  strSampleGroup="${strSampleGroup} GGH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} GGH_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} VBF_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2Nu2X_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo2L2Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HTo4Q_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_WWTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} ZH_HToLNuQQ_2LFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WminusH_HToWW_2LOSFilter_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Nu_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_ZZTo2L2Q_POWHEG"
  strSampleGroup="${strSampleGroup} WplusH_HToWW_2LOSFilter_POWHEG"
fi

declare -a SimSamples=( $(echo $strSampleGroup) )
for sample in "${SimSamples[@]}"; do
  for syst in "${systs[@]}"; do
    if [[ "$sample" != "GGH"* ]] && [[ "$syst" == "tHardJets"* ]]; then
      continue
    fi

    REQMEM=7168M
    JOBFLAV=tomorrow
    if [[ "$sample" == "ZH_HToLNuQQ_2LFilter_POWHEG" ]] || [[ "$sample" == "WminusH_ZZTo2L2Q_POWHEG" ]] || [[ "$sample" == "WplusH_ZZTo2L2Q_POWHEG" ]]; then
      REQMEM=2048M
      JOBFLAV=workday
    fi

    strargs="${arguments}"
    strargs="${strargs/<strSampleSet>/${sample}}"
    strargs="${strargs/<theGlobalSyst>/${syst}}"

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
  done
done
