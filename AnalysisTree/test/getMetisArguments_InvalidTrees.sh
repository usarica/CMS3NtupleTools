#!/bin/bash

indir=$1

rm -f tmpargslist.txt

declare -a arglist=( )
for f in $(grep -l -r -e "IFILE" ${indir}); do
  echo $f
  metis_dir=$(grep -e "Metis directory" ${f} | awk '{print $3}')
  metis_dir="../../NtupleMaker/test/tasks/${metis_dir}/logs/std_logs"
  for ifile in $(grep -e "IFILE" ${f} | awk '{print $2}'); do
    echo " - Checking ${ifile}..."
    tmploglist=( $(grep -l -r -e "IFILE: ${ifile}" ${metis_dir}) )
    for tmplf in "${tmploglist[@]}"; do
      echo "FOUND LOG FILE ${tmplf}"
    done
    tmpargs="$(grep -e 'args:' ${tmploglist[0]})"
    tmpargs="${tmpargs//'args: '}"
    echo ${tmpargs}
    arglist+=( "${tmpargs}" )
  done
done

for args in "${arglist[@]}"; do
  echo "\"${args}\" \\" >> tmpargslist.txt
done
