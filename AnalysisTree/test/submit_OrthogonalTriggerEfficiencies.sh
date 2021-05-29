#!/bin/bash

date=$1
period=$2
prodVersion=$3

useMETJERCorr=false
declare -i doSim=1
declare -i doData=1
for arg in "$@"; do
  if [[ "$arg" == "only_data" ]]; then
    doSim=0
    doData=1
  elif [[ "$arg" == "only_sim" ]]; then
    doSim=1
    doData=0
  elif [[ "$arg" == "useMETJERCorr="* ]]; then
    useMETJERCorr=${arg#*=}
  fi
done
if [[ "${useMETJERCorr}" != "true" ]] && [[ "${useMETJERCorr}" != "false" ]]; then
  echo "useMETJERCorr must be 'true' or 'false'"
fi

script=produceOrthogonalTriggerEfficiencies.cc
function=getTrees
jobdate="${date}_OrthogonalTriggerEfficiencies"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>",<ichunk>,<nchunks>,<theGlobalSyst>,<applyPUIdToAK4Jets>,<applyTightLeptonVetoIdToAK4Jets>,<use_MET_Puppi>,<use_MET_XYCorr>,<use_MET_JERCorr>,<use_MET_ParticleMomCorr>,<use_MET_p4Preservation>,<use_MET_corrections>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<prodVersion>/$prodVersion}"
arguments="${arguments/<applyPUIdToAK4Jets>/true}"
arguments="${arguments/<applyTightLeptonVetoIdToAK4Jets>/false}"
arguments="${arguments/<use_MET_Puppi>/false}"
arguments="${arguments/<use_MET_XYCorr>/true}"
arguments="${arguments/<use_MET_JERCorr>/$useMETJERCorr}"
arguments="${arguments/<use_MET_ParticleMomCorr>/true}"
arguments="${arguments/<use_MET_p4Preservation>/true}"
arguments="${arguments/<use_MET_corrections>/true}"

declare -a dataPeriods=( $period )
declare -a DataSampleList=( )
declare -a MCDataPeriods=( $period )
declare -a MCSysts=( sNominal )
declare -i dataYear

DataSampleList=( SingleLepton )

if [[ "$period" == "2018"* ]]; then
  dataYear=2018
  if [[ "$period" == "2018" ]]; then
    dataPeriods=( 2018A 2018B 2018C 2018D )
  fi
elif [[ "$period" == "2017"* ]]; then
  dataYear=2017
  if [[ "$period" == "2017" ]]; then
    dataPeriods=( 2017B 2017C 2017D 2017E 2017F )
  fi
elif [[ "$period" == "2016"* ]]; then
  dataYear=2016
  if [[ "$period" == "2016" ]]; then
    dataPeriods=( 2016B 2016C 2016D 2016E 2016F 2016G 2016H )
  fi
fi

MCDataPeriods=( $period )

csvfile="skimSamples_${dataYear}.csv"
maindir="/store/user/usarica/Offshell_2L2Nu"
maindir=$( ExecuteCompiledCommand GetStandardHostPathToStore ${maindir} t2.ucsd.edu )

if [[ $doSim -eq 1 ]]; then

for sample in $(readCMS3SkimSamplesFromCSV.py --csv=${csvfile} --sim --tree_req="SingleLepton"); do
  if [[ "$sample" != "/WJetsToLNu"* ]] || [[ "$sample" != *"amcatnlo"* ]]; then
    continue
  fi

  sampleDir=${sample//MINIAODSIM}

  echo "====="
  skimdir="${maindir}/Skims/${prodVersion}/${sampleDir}"
  proddir="${maindir}/Production/${prodVersion}/${sampleDir}"
  if [[ "$( ExecuteCompiledCommand DirectoryExists ${skimdir} )" == "false" ]]; then
    echo "$skimdir does not exist"
    skimdir="${maindir}/Skims/${prodVersion}_partial/${sampleDir}"
  fi
  if [[ "$( ExecuteCompiledCommand DirectoryExists ${skimdir} )" == "false" ]]; then
    echo "$skimdir does not exist"
    skimdir=$proddir
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
  if [[ $nchunks -lt 80 ]]; then
    nchunks=80
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
    sss=${spl}
    if [[ "${dataperiod}" == "2018"* ]] && [[ "${sss}" == "SingleElectron" ]]; then
       sss="EGamma"
    fi
    sample="${sss}_${dataperiod}"

    let nchunks=15
    REQMEM=4096M
    JOBFLAV=tomorrow
    if [[ "${dataperiod}" == "2017E" ]] || [[ "${dataperiod}" == "2017F" ]]; then
      let nchunks=30
    elif [[ "${dataperiod}" == "2018A" ]]; then
      let nchunks=30
    elif [[ "${dataperiod}" == "2018D" ]]; then
      let nchunks=50
    fi

    let ichunk=0
    while [[ $ichunk -lt $nchunks ]]; do
      strargs="${arguments}"
      strargs="${strargs/<strSampleSet>/${sample}}"
      strargs="${strargs/<ichunk>/${ichunk}}"
      strargs="${strargs/<nchunks>/${nchunks}}"
      strargs="${strargs/<theGlobalSyst>/sNominal}"
      strargs="${strargs/<period>/$dataperiod}"

      submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" memory="${REQMEM}" job_flavor="${JOBFLAV}"

      let ichunk=$ichunk+1
    done

  done
done

fi
