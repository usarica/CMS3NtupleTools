#!/bin/bash

date=$1
period=$2
prodVersion=$3
script=produceDileptonTriggerSurvey.cc
function=getTrees
jobdate="DileptonTriggerSurvey_${date}"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>",<ichunk>,<nchunks>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<prodVersion>/$prodVersion}"

declare -a MCDataPeriods=( $period )
declare -a MCSampleList=( )

DataSampleList=( Run )

if [[ "$period" == "2018"* ]]; then
  MCSampleList=( \
    /GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM \
    /VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM \
    /ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM \
    /TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM \
  )
elif [[ "$period" == "2017"* ]]; then
  MCSampleList=( \
    /GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM \
    /VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM \
    /ZZTo2L2Nu_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM \
    /TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM \
  )
elif [[ "$period" == "2016"* ]]; then
  MCSampleList=( \
    /GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM \
    /VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM \
    /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM \
    /ZZTo2L2Nu_13TeV_powheg_pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM \
    /TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM \
  )
fi

MCDataPeriods=( $period )


for sample in "${MCSampleList[@]}"; do
  sampleDir=${sample//MINIAODSIM}

  echo "====="
  skimdir="/ceph/cms/store/user/usarica/Offshell_2L2Nu/Skims/${prodVersion}/${sampleDir}"
  proddir="/ceph/cms/store/user/usarica/Offshell_2L2Nu/Production/${prodVersion}/${sampleDir}"
  skimdir=$( ExecuteCompiledCommand GetStandardHostPathToStore ${skimdir} t2.ucsd.edu )
  proddir=$( ExecuteCompiledCommand GetStandardHostPathToStore ${proddir} t2.ucsd.edu )
  if [[ "$( ExecuteCompiledCommand DirectoryExists ${skimdir} )" == "false" ]]; then
    echo "$skimdir does not exist"
    skimdir=$proddir
  fi
  if [[ "$( ExecuteCompiledCommand DirectoryExists ${skimdir} )" == "false" ]]; then
    continue
  fi

  let nfiles=$(ExecuteCompiledCommand lsdir $skimdir | grep -e ".root" | wc -l)
  echo "$sampleDir has $nfiles files"

  for dataperiod in "${MCDataPeriods[@]}"; do

    let ichunk=0
    while [[ $ichunk -lt $nfiles ]]; do

      strargs="${arguments}"
      strargs="${strargs/<strSampleSet>/${sample}}"
      strargs="${strargs/<ichunk>/${ichunk}}"
      strargs="${strargs/<nchunks>/${nfiles}}"
      strargs="${strargs/<period>/$dataperiod}"


      submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"

      let ichunk=$ichunk+1

    done
  done
  
  echo "====="
done
