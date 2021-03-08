#!/bin/bash

chkdir=$1
sleepdur=$2

while [[ 1 ]]; do
  for d in $(checkCMS3NtupleProduction.sh $chkdir | grep failed | awk '{print $1}' | sort); do
    resubmitCMS3NtupleProduction.sh $d
  done
  sleep $sleepdur
done
