#!/bin/bash

date=210504


#for d in DileptonEvents DileptonTriggerTnPEvents LGEvents LLGEvents SingleLeptonEvents SinglePhotonEGTnP SinglePhotonEvents 3LEvents LeptonTnP; do
for d in DileptonEvents LLGEvents SingleLeptonEvents SinglePhotonEvents 3LEvents; do
  scr="./submit_${d}.sh"

  for yy in 2016 2017 2018; do
    $scr ${date} $yy 201221_$yy all_imp_systs useMETJERCorr=false &
  done

  wait

  parallelizeCMS3Jobs.py --oldjobdir=output/${date}_${d}/*/* --newjobdir=output/${date}_${d}_parallel --nsubjobs=100 --nthreads=12
  resubmitCMS3NtupleProduction.sh output/${date}_${d}_parallel


  #for dd in $(checkCMS3NtupleProduction.sh output/${date}_${d} singleprod | grep -e failed | awk '{print $1}' | sort | uniq); do
  #  resubmitCMS3NtupleProduction.sh ${dd}
  #done
done
