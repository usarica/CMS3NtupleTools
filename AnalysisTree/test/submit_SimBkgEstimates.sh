#!/bin/bash

date=$1
period=$2
prodVersion=$3
script=produceSimBkgEstimates.cc
function=runDistributionsChain
jobdate="SimBkgEstimates_${date}"
arguments='"<period>","<prodVersion>","<strdate>",<theGlobalSyst>,<applyPUIdToAK4Jets>,<applyTightLeptonVetoIdToAK4Jets>,<use_MET_Puppi>,<use_MET_XYCorr>,<use_MET_JERCorr>,<use_MET_ParticleMomCorr>,<use_MET_p4Preservation>,<use_MET_corrections>'
arguments="${arguments/<period>/$period}"
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<prodVersion>/$prodVersion}"
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
    tEWDn tEWUp \
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
    eL1PrefiringDn eL1PrefiringUp \
    eTriggerEffDn eTriggerEffUp \
)


for syst in "${systs[@]}"; do
  strargs="${arguments}"
  strargs="${strargs/<theGlobalSyst>/${syst}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
done
