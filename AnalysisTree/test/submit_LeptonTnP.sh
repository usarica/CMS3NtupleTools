#!/bin/bash

date=$1
period=$2
prodVersion=$3
script=produceLeptonEfficiencies.cc
function=getTrees
jobdate="LeptonTnP_${date}"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>",<ichunk>,<nchunks>,<theGlobalSyst>,<applyPUIdToAK4Jets>,<applyTightLeptonVetoIdToAK4Jets>,<vetoExtraNonOverlappingLeptons>,<hardProcessFallback>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<prodVersion>/$prodVersion}"
arguments="${arguments/<applyPUIdToAK4Jets>/true}"
arguments="${arguments/<applyTightLeptonVetoIdToAK4Jets>/false}"

declare -a csvList

declare -a dataPeriods
declare -a MCDataPeriods
declare -a DataSampleList
declare -a MCSampleList

DataSampleList=( EGamma SingleMuon )

if [[ "$period" == "2018" ]]; then
  dataPeriods=( 2018A 2018B 2018C 2018D )
  MCSampleList=( \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM \
    /ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM \
    /ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM \
    /ZZTo2L2Nu_mZMin-18_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM \
    /ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM \
  )
elif [[ "$period" == "2017" ]]; then
  dataPeriods=( 2017B 2017C 2017D 2017E 2017F )
  MCSampleList=( \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM \
    /ZZTo2L2Nu_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM \
    /ZZTo2L2Nu_mZMin-18_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM \
    /ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM \
    /ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM \
  )
elif [[ "$period" == "2016" ]]; then
  dataPeriods=( 2016B 2016C 2016D 2016E 2016F 2016G 2016H )
  MCSampleList=( \
    /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM \
    /ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM \
    /ZZTo2L2Nu_13TeV_powheg_pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM \
    /ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM \
    /ZZTo4L_13TeV_powheg_pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM \
  )
fi

#MCDataPeriods=("${dataPeriods[@]}")
MCDataPeriods=( $period )


for sample in "${MCSampleList[@]}"; do
  sampleDir=${sample//MINIAODSIM}

  echo "====="
  skimdir="/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/${prodVersion}/${sampleDir}"
  skimdir=$( ExecuteCompiledCommand GetStandardHostPathToStore ${skimdir} t2.ucsd.edu )
  if [[ "$( ExecuteCompiledCommand DirectoryExists ${skimdir} )" == "false" ]]; then
    echo "$skimdir does not exist"
    continue
  fi

  let nfiles=$(ExecuteCompiledCommand lsdir $skimdir | grep -e ".root" | wc -l)
  echo "$sampleDir has $nfiles files"

  for dataperiod in "${MCDataPeriods[@]}"; do
    for syst in sNominal ePUDn ePUUp; do

      let ichunk=0
      while [[ $ichunk -lt $nfiles ]]; do
        strargs="${arguments}"
        strargs="${strargs/<strSampleSet>/${sample}}"
        strargs="${strargs/<ichunk>/${ichunk}}"
        strargs="${strargs/<nchunks>/${nfiles}}"
        strargs="${strargs/<theGlobalSyst>/${syst}}"
        strargs="${strargs/<period>/$dataperiod}"
        if [[ "${sample}" == "/ZZTo2L2Nu"* ]];then
          strargs="${strargs/<vetoExtraNonOverlappingLeptons>/true}"
          #strargs="${strargs/<hardProcessFallback>/true}"
          strargs="${strargs/<hardProcessFallback>/false}"
        elif [[ "${sample}" == "/ZZTo4L"* ]];then
          strargs="${strargs/<vetoExtraNonOverlappingLeptons>/false}"
          #strargs="${strargs/<hardProcessFallback>/true}"
          strargs="${strargs/<hardProcessFallback>/false}"
        else
          strargs="${strargs/<vetoExtraNonOverlappingLeptons>/true}"
          strargs="${strargs/<hardProcessFallback>/false}"
        fi

        submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
        let ichunk=$ichunk+1
      done

    done
  done
  
  echo "====="
done

for spl in "${DataSampleList[@]}"; do
  for dataperiod in "${dataPeriods[@]}"; do
    sample="${spl}_${dataperiod}"

    strargs="${arguments}"
    strargs="${strargs/<strSampleSet>/${sample}}"
    strargs="${strargs/<ichunk>/0}"
    strargs="${strargs/<nchunks>/0}"
    strargs="${strargs/<theGlobalSyst>/sNominal}"
    strargs="${strargs/<period>/$dataperiod}"
    strargs="${strargs/<vetoExtraNonOverlappingLeptons>/true}"
    strargs="${strargs/<hardProcessFallback>/false}"

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
  done
done
