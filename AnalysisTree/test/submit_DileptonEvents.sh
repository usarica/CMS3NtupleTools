#!/bin/bash

date=$1
period=$2
prodVersion=$3
script=produceDileptonEvents.cc
function=getTrees
jobdate="DileptonEvents_${date}"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>",<ichunk>,<nchunks>,<theGlobalSyst>,<computeMEs>,<applyPUIdToAK4Jets>,<applyTightLeptonVetoIdToAK4Jets>,<use_MET_Puppi>,<use_MET_XYCorr>,<use_MET_JERCorr>,<use_MET_ParticleMomCorr>,<use_MET_p4Preservation>,<use_MET_corrections>'
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<prodVersion>/$prodVersion}"
arguments="${arguments/<computeMEs>/true}"
arguments="${arguments/<applyPUIdToAK4Jets>/true}"
arguments="${arguments/<applyTightLeptonVetoIdToAK4Jets>/false}"
arguments="${arguments/<use_MET_Puppi>/false}"
arguments="${arguments/<use_MET_XYCorr>/true}"
arguments="${arguments/<use_MET_JERCorr>/true}"
arguments="${arguments/<use_MET_ParticleMomCorr>/true}"
arguments="${arguments/<use_MET_p4Preservation>/true}"
arguments="${arguments/<use_MET_corrections>/true}"

declare -a dataPeriods=( $period )
declare -a DataSampleList=( )
declare -a MCDataPeriods=( $period )
declare -a MCSysts=( sNominal )
declare -i dataYear

DataSampleList=( Run )

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
for sample in $(readCMS3SkimSamplesFromCSV.py --csv=${csvfile} --sim --tree_req="Dilepton"); do
  sampleDir=${sample//MINIAODSIM}

  echo "====="
  skimdir="/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/${prodVersion}/${sampleDir}"
  proddir="/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Production/${prodVersion}/${sampleDir}"
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

for spl in "${DataSampleList[@]}"; do
  for dataperiod in "${dataPeriods[@]}"; do
    sample="${spl}${dataperiod}"

    let nchunks=4
    REQMEM=4096M
    JOBFLAV=tomorrow
    if [[ "${dataperiod}" == "2017E" ]] || [[ "${dataperiod}" == "2017F" ]]; then
      JOBFLAV=testmatch
      let nchunks=8
    elif [[ "${dataperiod}" == "2018A" ]]; then
      JOBFLAV=testmatch
      let nchunks=10
    elif [[ "${dataperiod}" == "2018D" ]]; then
      JOBFLAV=testmatch
      let nchunks=30
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
