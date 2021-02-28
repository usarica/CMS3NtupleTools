#!/bin/bash

date=$1
period=$2
prodVersion=$3

useMETJERCorr=false
declare -i doSim=1
declare -i doData=1
declare -i doSysts=0 # Do systematics
for arg in "$@"; do
  if [[ "$arg" == "only_data" ]]; then
    doSim=0
    doData=1
  elif [[ "$arg" == "only_sim" ]]; then
    doSim=1
    doData=0
  elif [[ "$arg" == "all_systs" ]] || [[ "$arg" == "all_imp_systs" ]]; then
    doSysts=1
  elif [[ "$arg" == "useMETJERCorr="* ]]; then
    useMETJERCorr=${arg#*=}
  fi
done
if [[ "${useMETJERCorr}" != "true" ]] && [[ "${useMETJERCorr}" != "false" ]]; then
  echo "useMETJERCorr must be 'true' or 'false'"
fi

script=produceSinglePhotonEGTnP.cc
function=getTrees
jobdate="${date}_SinglePhotonEGTnP"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>",<ichunk>,<nchunks>,<theGlobalSyst>,<hardProcessFallback>,<applyPUIdToAK4Jets>,<applyTightLeptonVetoIdToAK4Jets>,<use_MET_Puppi>,<use_MET_XYCorr>,<use_MET_JERCorr>,<use_MET_ParticleMomCorr>,<use_MET_p4Preservation>,<use_MET_corrections>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<prodVersion>/$prodVersion}"
arguments="${arguments/<hardProcessFallback>/false}"
arguments="${arguments/<applyPUIdToAK4Jets>/true}"
arguments="${arguments/<applyTightLeptonVetoIdToAK4Jets>/false}"
arguments="${arguments/<use_MET_Puppi>/false}"
arguments="${arguments/<use_MET_XYCorr>/true}"
arguments="${arguments/<use_MET_JERCorr>/$useMETJERCorr}"
arguments="${arguments/<use_MET_ParticleMomCorr>/true}"
arguments="${arguments/<use_MET_p4Preservation>/true}"
arguments="${arguments/<use_MET_corrections>/true}"

declare -a dataPeriods=( )
declare -a MCDataPeriods=( $period )
declare -a DataSampleList=( EGamma )
declare -a MCSampleList=( )

if [[ "$period" == "2018" ]]; then
  dataPeriods=( 2018A 2018B 2018C 2018D )
  MCSampleList=( \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM \
  )
elif [[ "$period" == "2017" ]]; then
  dataPeriods=( 2017B 2017C 2017D 2017E 2017F )
  MCSampleList=( \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM \
    /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM \
  )
elif [[ "$period" == "2016" ]]; then
  dataPeriods=( 2016B 2016C 2016D 2016E 2016F 2016G 2016H )
  MCSampleList=( \
    /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM \
  )
fi

declare -a MCSysts=( sNominal )
# Does not matter which flag is used...
if [[ $doSysts -eq 1 ]]; then
  MCSysts+=( \
    ePUDn ePUUp \
  )
fi


maindir="/store/user/usarica/Offshell_2L2Nu"
maindir=$( ExecuteCompiledCommand GetStandardHostPathToStore ${maindir} t2.ucsd.edu )

if [[ $doSim -eq 1 ]]; then

for sample in "${MCSampleList[@]}"; do
  sampleDir=${sample//MINIAODSIM}

  echo "====="
  skimdir="${maindir}/Skims/${prodVersion}/${sampleDir}"
  if [[ "$( ExecuteCompiledCommand DirectoryExists ${skimdir} )" == "false" ]]; then
    echo "$skimdir does not exist"
    skimdir="${maindir}/Skims/${prodVersion}_partial/${sampleDir}"
  fi
  if [[ "$( ExecuteCompiledCommand DirectoryExists ${skimdir} )" == "false" ]]; then
    continue
  fi

  let nfiles=$(ExecuteCompiledCommand lsdir $skimdir | grep -e ".root" | wc -l)
  echo "$sampleDir has $nfiles files"

  nchunks=$((nfiles / 4))
  nfeff=$((nchunks * 4))
  if [[ $nfeff -lt $nfiles ]];then
    nchunks=$((nchunks + 1))
  fi
  echo " - There will be $nchunks chunks"
  if [[ $nchunks -eq 0 ]]; then
    echo "Oh-uh..."
  fi

  for dataperiod in "${MCDataPeriods[@]}"; do
    for syst in "${MCSysts[@]}"; do

      let ichunk=0
      while [[ $ichunk -lt $nchunks ]]; do
        strargs="${arguments}"
        strargs="${strargs/<strSampleSet>/${sample}}"
        strargs="${strargs/<ichunk>/${ichunk}}"
        strargs="${strargs/<nchunks>/${nchunks}}"
        strargs="${strargs/<theGlobalSyst>/${syst}}"
        strargs="${strargs/<period>/$dataperiod}"

        submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
        let ichunk=$ichunk+1
      done

    done
  done
  
  echo "====="
done

fi

if [[ $doData -eq 1 ]]; then

for spl in "${DataSampleList[@]}"; do
  for dataperiod in "${dataPeriods[@]}"; do
    sample="${spl}_${dataperiod}"

    strargs="${arguments}"
    strargs="${strargs/<strSampleSet>/${sample}}"
    strargs="${strargs/<ichunk>/0}"
    strargs="${strargs/<nchunks>/0}"
    strargs="${strargs/<theGlobalSyst>/sNominal}"
    strargs="${strargs/<period>/$dataperiod}"

    REQMEM=2048M
    JOBFLAV=tomorrow
    if [[ "${dataperiod}" == "2017E" ]] || [[ "${dataperiod}" == "2017F" ]] || [[ "${dataperiod}" == "2018A" ]] || [[ "${dataperiod}" == "2018D" ]]; then
      REQMEM=4096M
      JOBFLAV=testmatch
    fi

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"
  done
done

fi
