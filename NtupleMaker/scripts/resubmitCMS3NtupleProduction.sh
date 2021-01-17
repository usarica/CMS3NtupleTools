#!/bin/bash

chkdir=$1
declare -i jobprio=0
for farg in "$@"; do
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "${fargl}" == "priority="* ]]; then
    jobprio=${fargl//'priority='}
  fi
done


for f in $(find $chkdir -name condor.sub); do
  d=${f//\/condor.sub}
  cd $d
  echo "Processing $d"
  rm -f Logs/prior_record.tar

  for prevjob in $(ls ./ | grep ".log"); do
    prevjob=${prevjob//".log"}
    strstdlogs=""
    for stdlog in $(ls Logs | grep -e ${prevjob}); do
      strstdlogs="${strstdlogs} Logs/${stdlog}"
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
done
