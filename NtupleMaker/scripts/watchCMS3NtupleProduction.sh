#!/bin/bash

declare -r chkdir=$1
declare -i sleepdur=1800
check_opts=""
resubmit_opts=""
for farg in "$@"; do
  fargl="$(echo $farg | awk '{print tolower($0)}')"

  if [[ "$fargl" == "sleep="* ]]; then
    let sleepdur=${fargl//'sleep='}
  elif [[ "$fargl" == "check-opts="* ]]; then
    check_opts="${farg#*=}"
  elif [[ "$fargl" == "resubmit-opts="* ]]; then
    resubmit_opts="${farg#*=}"
  fi
done

while [[ 1 ]]; do
  for d in $(checkCMS3NtupleProduction.sh ${chkdir} ${check_opts} | grep failed | awk '{print $1}' | sort); do
    resubmitCMS3NtupleProduction.sh ${d} ${resubmit_opts}
  done
  sleep ${sleepdur}
done
