#!/bin/bash

declare -r chkdir="$1"
declare -i jobprio=0
declare -i maxparallel=1
for farg in "$@"; do
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "${fargl}" == "priority="* ]]; then
    jobprio=${fargl//'priority='}
  elif [[ "$fargl" == "maxparallel="* ]]; then
    let maxparallel=${fargl//'maxparallel='}
  fi
done


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
resubmitDirectory(){
  local f=$1
  local d=${f//\/condor.sub}
  cd $d
  echo "Processing $d"
  rm -f Logs/prior_record.tar

  for prevjob in $(ls ./ | grep ".log"); do
    local prevjob=${prevjob//".log"}
    strstdlogs=""
    for stdlog in $(ls Logs | grep -e ${prevjob}); do
      local strstdlogs="${strstdlogs} Logs/${stdlog}"
    done
    tar Jcf "prior_record.${prevjob}.tar" ${strstdlogs} "${prevjob}.log" --exclude={*.tar}
    rm "${prevjob}.log"
    for stdlog in $(echo ${strstdlogs}); do
      rm ${stdlog}
    done
  done

  condor_submit condor.sub

  if [[ $jobprio -ne 0 ]]; then
    for prevjob in $(ls ./ | grep ".log"); do
      prevjob=${prevjob//".log"}
      condor_prio $prevjob -p $jobprio
    done
  fi

  cd - &> /dev/null
}


for f in $(find $chkdir -name condor.sub); do
  if [[ $maxparallel -le 1 ]]; then
    resubmitDirectory $f
  else
    resubmitDirectory $f &
    job_limit $maxparallel
  fi
done
