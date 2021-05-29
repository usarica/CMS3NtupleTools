#!/bin/bash

date=$1
period=$2
prodVersion=$3
ntupleVersion=$4

script=produceNRBEstimates.cc
function=runDistributionsChain
jobdate="${date}_NRBEstimates"
arguments='"<period>","<prodVersion>","<ntupleVersion>","<strdate>",<theGlobalSyst>,<fast_mode>,<applyPUIdToAK4Jets>,<applyTightLeptonVetoIdToAK4Jets>,<use_MET_Puppi>,<use_MET_XYCorr>,<use_MET_JERCorr>,<use_MET_ParticleMomCorr>,<use_MET_p4Preservation>,<use_MET_corrections>'
arguments="${arguments/<period>/$period}"
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<prodVersion>/$prodVersion}"
arguments="${arguments/<ntupleVersion>/$ntupleVersion}"
arguments="${arguments/<fast_mode>/false}"
arguments="${arguments/<applyPUIdToAK4Jets>/true}"
arguments="${arguments/<applyTightLeptonVetoIdToAK4Jets>/false}"
arguments="${arguments/<use_MET_Puppi>/false}"
arguments="${arguments/<use_MET_XYCorr>/true}"
arguments="${arguments/<use_MET_JERCorr>/false}"
arguments="${arguments/<use_MET_ParticleMomCorr>/true}"
arguments="${arguments/<use_MET_p4Preservation>/true}"
arguments="${arguments/<use_MET_corrections>/true}"

# Run only on those that change data and MC alike
declare -a systs=( \
  sNominal \
  eEleEffStatDn eEleEffStatUp \
  eEleEffSystDn eEleEffSystUp \
  eMuEffStatDn eMuEffStatUp \
  eMuEffSystDn eMuEffSystUp \
  eTriggerEffDn eTriggerEffUp \
)


for syst in "${systs[@]}"; do
  strargs="${arguments}"
  strargs="${strargs/<theGlobalSyst>/${syst}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
done
