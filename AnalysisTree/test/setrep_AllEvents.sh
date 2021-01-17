#!/bin/bash

date=210107


while [[ 1 ]]; do
  for d in DileptonEvents DileptonTriggerTnPEvents LGEvents LLGEvents SingleLeptonEvents SinglePhotonEGTnP SinglePhotonEvents 3LEvents LeptonEfficiencies; do
    hadoop fs -setrep -R 4 /cms/store/user/usarica/Offshell_2L2Nu/Worker/output/${d}/SkimTrees/${date}
  done
  sleep 7200
done
