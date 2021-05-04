#!/bin/bash

job_limit(){
  local joblist=( )
  # Test for single positive integer input
  if [[ $# -eq 1 ]] && [[ $1 =~ ^[1-9][0-9]*$ ]]
  then
    # Check number of running jobs
    joblist=( $(jobs -rp) )
    while [[ ${#joblist[*]} -ge $1 ]]; do
      # Wait for any job to finish
      command='wait '${joblist[0]}
      for job in ${joblist[@]:1}; do
        command+=' || wait '$job
      done
      eval ${command}
      joblist=( $(jobs -rp) )
    done
  fi
}

date=$1
period=$2
fitVersion=$3
nthreads=1
if [[ "$4" != "" ]]; then
  nthreads=$4
fi

script=produceLeptonEfficiencies_RooFit.cc
function=getNLLRecovery
jobdate="${date}_LeptonTnP_FitNLL_${period}"
arguments='"<strdate>","<period>","<fitVersion>","<wsdcname>"'
arguments="${arguments/<period>/$period}"
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<fitVersion>/$fitVersion}"

maindir="/store/user/usarica/Offshell_2L2Nu/Worker/output/LeptonEfficiencies/DataFits/${fitVersion}/${period}"
maindir=$( ExecuteCompiledCommand GetStandardHostPathToStore ${maindir} t2.ucsd.edu )

wsdcnames=( $(ExecuteCompiledCommand lsdir ${maindir}) )
for wsdcname in "${wsdcnames[@]}"; do
  skimdir=${maindir}/${wsdcname}
  if [[ "$( ExecuteCompiledCommand FileExists ${skimdir}/combined_withSnapshot.root )" == "false" ]] && [[ "$( ExecuteCompiledCommand FileExists ${skimdir}/combined_withSnapshot_withRobustFit.root )" == "false" ]]; then
    continue
  fi

  strargs="${arguments}"
  strargs="${strargs/<wsdcname>/${wsdcname}}"

  submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}" &
  job_limit $nthreads
done

wait
