#!/bin/bash

date=210224


while [[ 1 ]]; do
#for d in DileptonEvents DileptonTriggerTnPEvents LGEvents LLGEvents SingleLeptonEvents SinglePhotonEGTnP SinglePhotonEvents 3LEvents LeptonTnP; do
for d in DileptonEvents LGEvents LLGEvents SingleLeptonEvents SinglePhotonEvents 3LEvents; do
  scr="./submit_${d}.sh"

  #$scr ${date} 2018 201221_2018 only_sim useMETJERCorr=false
  #$scr ${date} 2017 201221_2017 only_sim useMETJERCorr=false
  #$scr ${date} 2016 201221_2016 only_sim useMETJERCorr=false

  for dd in $(checkCMS3NtupleProduction.sh output/${date}_${d} singleprod | grep -e failed | awk '{print $1}' | sort | uniq); do
    resubmitCMS3NtupleProduction.sh ${dd}
  done
done
sleep 1800
done
