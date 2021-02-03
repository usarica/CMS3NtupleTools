#!/bin/bash

chkdir=$1
let multiprod=1
let islongfile=0
let skiplongfile=0
if [[ "$2" != "" ]];then
  if [[ "$2" == "singleprod" ]]; then
    let multiprod=0
  elif [[ "$2" == "longfile" ]]; then
    let islongfile=1
  elif [[ "$2" == "skiplongfile" ]]; then
    let skiplongfile=1
  fi
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

declare -a runningjobs=( $(condor_q -af:j '') )

for f in $(find $chkdir -name condor.sub); do
  d=${f//\/condor.sub}
  if [[ ! -d $d ]];then
    continue
  fi

  let countOK=0
  let dirok=1
  let nsubjobs=0
  let nRunningJobs=0
  for joblog in $(ls $d | grep -e ".log"); do
    jobnumber=${joblog//".log"}
    logfilename="log_job.${jobnumber}.txt"

    let nsubjobs=$nsubjobs+1

    resb=( )
    rese=( )
    resf=( )
    ress=( )

    let job_is_running=0
    for jobid in "${runningjobs[@]}"; do
      if [[ "${jobid}" == "${jobnumber}" ]]; then
        let job_is_running=1
        break
      fi
    done
    
    if [[ ${job_is_running} -eq 1 ]]; then
      echo "Job $jobnumber for $d is still running"
      let nUNKNOWN=$nUNKNOWN+1
      let dirok=0
      let nRunningJobs=${nRunningJobs}+1
      continue
    fi

    fread="$d/Logs/$logfilename"

    if [[ -e ${fread} ]]; then
      if [[ $skiplongfile -eq 1 ]]; then
        let freadsize=$(stat --format=%s $fread)
        if [[ $freadsize -gt 1000000 ]]; then
          echo "Skipping $fread because its size = $freadsize > 1000000"
          continue
        fi
        #echo "File $fread with size $freadsize is not long..."
      fi
    fi

    if [[ -e ${fread} ]]; then
      if [[ $islongfile -eq 1 ]]; then
        tmpfile="tail_$logfilename"
        echo "Truncating $fread to $tmpfile"
        tail -1000 $fread &> $tmpfile
        fread=$tmpfile
        echo "Will read $fread"
      fi
    fi

    if [[ -e ${fread} ]]; then
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
      done < "$fread"
    fi

    if [[ $islongfile -eq 1 ]]; then
      rm -f $fread
    fi

    let size_resb=${#resb[@]}
    let size_rese=${#rese[@]}
    let size_resf=${#resf[@]}
    let size_ress=${#ress[@]}

    if [[ $size_resb -gt 0 ]] && [[ $size_resb -eq $size_rese ]] && [[ $size_resb -eq $size_ress ]] && [[ $size_resb -eq $size_resf ]];then
      let nOutputExist=0
      for rf in "${resf[@]}";do
        rf="${rf//*'gsiftp://gftp.t2.ucsd.edu'}"

        if [[ -s $rf ]];then
          let nOutputExist=${nOutputExist}+1
        fi
      done
      if [[ $nOutputExist -eq $size_resf ]];then
        echo "Job $jobnumber for $d ran successfully with $nOutputExist files."
        let nOK=$nOK+1
        let countOK=$countOK+1
      else
        echo "Job $jobnumber for $d ran successfully, but some files do not exist. (Nbegin, Ncopyrun, Nsuccess, Nend, Nexists) = ( $size_resb, $size_resf, $size_ress, $size_rese, $nOutputExist )"
        let nFILEDNE=$nFILEDNE+1
        let dirok=0
      fi
    else
      echo "Job $jobnumber for $d did not run successfully. (Nbegin, Ncopyrun, Nsuccess, Nend) = ( $size_resb, $size_resf, $size_ress, $size_rese )"
      let nFAIL=$nFAIL+1
      let dirok=0
    fi

  done

  if [[ $countOK -gt 0 ]] && [[ $dirok -eq 0 ]];then
    if [[ $multiprod -eq 1 ]];then
      echo "$d has multiple submissions with $countOK / $nsubjobs success rate, but all subjobs were required to succeed."
    else
      echo "$d has multiple submissions with $countOK / $nsubjobs success rate, but the user specified the folders to be treated as single jobs."
      let dirok=1
    fi
  fi

  if [[ $nsubjobs -eq 0 ]];then
    echo "$d does not have any subjobs run yet."
    let nFAIL=$nFAIL+1
    let dirok=0
  fi

  if [[ $dirok -eq 1 ]];then
    TARFILE="${d}.tar"
    TARFILEFirst=${TARFILE%/*}
    TARFILELast=${TARFILE##*/}
    if [[ ${#TARFILELast} -gt 255 ]]; then
      TARFILE=${TARFILEFirst}/${TARFILELast//_}
    fi
    rm -f $TARFILE
    tar Jcf ${TARFILE} $d --exclude={*.tar}
    if [[ $? -eq 0 ]];then
      echo "- Compressed successfully, so removing the directory"
      rm -rf $d
    fi
  elif [[ $nRunningJobs -gt 0 ]]; then
    echo "$d is still running."
  else
    echo "$d failed."
  fi

done

echo "(OK:COPY_FAIL:FILE_DNE:FAIL:UNKNOWN) = (${nOK}:${nCOPYFAIL}:${nFILEDNE}:${nFAIL}:${nUNKNOWN})"
