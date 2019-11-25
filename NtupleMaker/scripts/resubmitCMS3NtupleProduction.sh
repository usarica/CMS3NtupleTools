#!/bin/bash

chkdir=$1

for f in $(find $chkdir -name condor.sub); do
  d=${f//\/condor.sub}
  cd $d
  echo "Processing $d"
  rm -f Logs/prior_record.tar
  prevjob=$(ls ./ | grep ".log")
  if [[ ! -z $prevjob ]];then
    prevjob=${prevjob//".log"}
    tar Jcf "prior_record.${prevjob}.tar" Logs/* "${prevjob}.log" --exclude={*.tar}
    rm "${prevjob}.log"
    rm Logs/*
  else
    echo "Could not identify the prior job number in folder $d"
  fi

  condor_submit condor.sub

  cd - &> /dev/null
done
