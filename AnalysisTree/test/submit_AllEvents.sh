#!/bin/bash

date=210107


while [[ 1 ]]; do
for d in DileptonEvents DileptonTriggerTnPEvents LGEvents LLGEvents SingleLeptonEvents SinglePhotonEGTnP SinglePhotonEvents 3LEvents LeptonTnP; do
  scr="./submit_${d}.sh"

  #$scr ${date} 2018 200906_2018 all_imp_systs only_sim useMETJERCorr=false
  #$scr ${date} 2017 201221_2017 all_imp_systs only_sim useMETJERCorr=false
  #$scr ${date} 2016 200906_2016 all_imp_systs only_sim useMETJERCorr=false

  for dd in $(checkCMS3NtupleProduction.sh output/${date}_${d} | grep -e failed | awk '{print $1}' | sort | uniq); do
    resubmitCMS3NtupleProduction.sh ${dd}
  done
done
sleep 3600
done
