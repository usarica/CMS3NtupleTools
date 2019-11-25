#!/bin/bash

chkdir=$1
let multiprod=1
if [[ "$2" != "" ]];then
  let multiprod=$2
fi


let nOK=0
let nCOPYFAIL=0
let nFAIL=0
let nFILEDNE=0
let nUNKNOWN=0

declare -a resb
declare -a rese
declare -a resf
declare -a ress

cd $chkdir

for d in $(ls ./); do
  if [[ ! -d $d ]];then
    continue
  fi

  let countOK=0
  let dirok=1
  let nsubjobs=0
  for logfilename in $(ls ./$d/Logs | grep -e "log_"); do
    let nsubjobs=$nsubjobs+1

    resb=( )
    rese=( )
    resf=( )
    ress=( )

    while IFS='' read -r line || [[ -n "$line" ]]; do
      if [[ "$line" == *"begin copying output"* ]]; then
        resb+=( "$line" )
      elif [[ "$line" == *"end copying output"* ]]; then
        rese+=( "$line" )
      elif [[ "$line" == *"Running: env -i "* ]]; then
        resf+=( "$line" )
      elif [[ "$line" == *"Copied successfully"* ]]; then
        ress+=( "$line" )
      fi
    done < "$d/Logs/$logfilename"

    let size_resb=${#resb[@]}
    let size_rese=${#rese[@]}
    let size_resf=${#resf[@]}
    let size_ress=${#ress[@]}

    if [[ ! -s $d/Logs/$logfilename ]];then
      echo "$d status is not yet determined"
      let nUNKNOWN=$nUNKNOWN+1
      let dirok=0
    elif [[ $size_resb -gt 0 ]] && [[ $size_resb -eq $size_rese ]] && [[ $size_resb -eq $size_ress ]] && [[ $size_resb -eq $size_resf ]];then
      let nOutputExist=0
      for rf in "${resf[@]}";do
        rf="${rf//*'gsiftp://gftp.t2.ucsd.edu'}"

        if [[ -s $rf ]];then
          let nOutputExist=${nOutputExist}+1
        fi
      done
      if [[ $nOutputExist -eq $size_resf ]];then
        echo "$d ran successfully with $nOutputExist files."
        let nOK=$nOK+1
        let countOK=$countOK+1
      else
        echo "$d ran successfully, but some files do not exist! (Nbegin, Ncopyrun, Nsuccess, Nend, Nexists) = ( $size_resb, $size_resf, $size_ress, $size_rese, $nOutputExist )"
        let nFILEDNE=$nFILEDNE+1
        let dirok=0
      fi
    else
      echo "$d failed. (Nbegin, Ncopyrun, Nsuccess, Nend) = ( $size_resb, $size_resf, $size_ress, $size_rese )"
      let nFAIL=$nFAIL+1
      let dirok=0
    fi

  done

  if [[ $countOK -gt 0 ]] && [[ $dirok -eq 0 ]];then
    if [[ $multiprod -eq 1 ]];then
      echo "$d has multiple submissions with $countOK / $nsubjobs success rate, but the folder will be treated as if it failed."
    else
      echo "$d has multiple submissions with $countOK / $nsubjobs success rate, but the user specified the folders to be treated as single jobs."
      let dirok=1
    fi
  fi

  if [[ $dirok -eq 1 ]];then
    TARFILE="${d}.tar"
    rm -f $TARFILE
    tar Jcf ${TARFILE} $d --exclude={*.tar}
    if [[ $? -eq 0 ]];then
      echo "- Compressed successfully, so removing the directory"
      rm -rf $d
    fi
  fi

done

cd -

echo "(OK:COPY_FAIL:FILE_DNE:FAIL:UNKNOWN) = (${nOK}:${nCOPYFAIL}:${nFILEDNE}:${nFAIL}:${nUNKNOWN})"
