#!/bin/bash

chkdir=$1
sleepdur=1800
if [[ "$2" != "" ]]; then
  sleepdur=$2
fi

while [[ 1 ]]; do
  for d in $(checkCMS3NtupleProduction.sh $chkdir | grep failed | awk '{print $1}' | sort); do
    resubmitCMS3NtupleProduction.sh $d
  done
  sleep $sleepdur
done
