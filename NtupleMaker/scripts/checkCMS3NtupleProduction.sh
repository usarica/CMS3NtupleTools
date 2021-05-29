#!/bin/bash

declare -r chkdir="$1"
declare -i multiprod=1
declare -i islongfile=0
declare -i skiplongfile=0
declare -i skipcondorcheck=0
declare -i fast_tar=0
declare -i no_tar=0
declare -i maxparallel=1
for farg in "$@"; do
  fargl="$(echo $farg | awk '{print tolower($0)}')"

  if [[ "$fargl" == "singleprod" ]]; then
    let multiprod=0
  elif [[ "$fargl" == "longfile" ]]; then
    let islongfile=1
  elif [[ "$fargl" == "skiplongfile" ]]; then
    let skiplongfile=1
  elif [[ "$fargl" == "skipcondorcheck" ]]; then
    let skipcondorcheck=1
  elif [[ "$fargl" == "fast-tar" ]]; then
    let fast_tar=1
  elif [[ "$fargl" == "no-tar" ]]; then
    let no_tar=1
  elif [[ "$fargl" == "maxparallel="* ]]; then
    let maxparallel=${fargl//'maxparallel='}
  fi
done


declare -i nOK=0
declare -i nCOPYFAIL=0
declare -i nFAIL=0
declare -i nFILEDNE=0
declare -i nUNKNOWN=0

declare -a resb
declare -a rese
declare -a resf
declare -a ress

declare -a runningjobs=( )
if [[ $skipcondorcheck -eq 0 ]]; then
  runningjobs=( $(condor_q -c 'JobStatus<=2 || JobStatus>=5' -af:j '') )
fi

job_limit(){
  local joblist=( )
  # Test for single positive integer input
  if [[ $# -eq 1 ]] && [[ $1 =~ ^[1-9][0-9]*$ ]]
  then
    # Check number of running jobs
    joblist=( $(jobs -rp) )
    while [[ ${#joblist[*]} -ge $1 ]]; do
      # Wait for any job to finish
      command='wait '${joblist[0]}
      for job in ${joblist[@]:1}; do
        command+=' || wait '$job
      done
      eval ${command}
      joblist=( $(jobs -rp) )
    done
  fi
}

tarAndRemoveDirectory(){
  local d=$1
  if [[ ${no_tar} -eq 0 ]]; then
    local TARFILE="${d}.tar"
    local TARFILEFirst=${TARFILE%/*}
    local TARFILELast=${TARFILE##*/}
    if [[ ${#TARFILELast} -gt 255 ]]; then
      TARFILE=${TARFILEFirst}/${TARFILELast//_}
    fi
    rm -f ${TARFILE}
    if [[ ${fast_tar} -eq 0 ]]; then
      tar Jcf ${TARFILE} $d --exclude={*.tar}
    else
      tar cf ${TARFILE} $d --exclude={*.tar}
    fi
    if [[ $? -eq 0 ]]; then
      echo "- Compressed successfully, so removing the directory"
      rm -rf $d
    fi
  fi
}

checkDirectory(){
  local f=$1
  local d=${f//\/condor.sub}
  if [[ ! -d $d ]];then
    continue
  fi

  local countOK=0
  local dirok=1
  local nsubjobs=0
  local nRunningJobs=0
  local job_is_running=0
  local size_resb=0
  local size_rese=0
  local size_resf=0
  local size_ress=0
  for joblog in $(ls $d | grep -e ".log"); do
    local jobnumber=${joblog//".log"}
    local logfilename="log_job.${jobnumber}.txt"

    let nsubjobs=$nsubjobs+1

    local resb=( )
    local rese=( )
    local resf=( )
    local ress=( )

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

    local fread="$d/Logs/$logfilename"

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
        local tmpfile="$d/Logs/tail_$logfilename"
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
        elif [[ "$line" == *"OUTPUTFILE: "* ]]; then
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

    if [[ $size_resb -gt 0 ]] && [[ $size_resb -eq $size_rese ]] && [[ $size_resb -eq $size_ress ]] && [[ $size_resb -eq $size_resf ]]; then
      {
      local nOutputExist=0
      for rf in "${resf[@]}";do
        rf="${rf//*'OUTPUTFILE: '}"

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
      }
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
    tarAndRemoveDirectory $d
  elif [[ $nRunningJobs -gt 0 ]]; then
    echo "$d is still running."
  else
    echo "$d failed."
  fi
}

for f in $(find $chkdir -name condor.sub); do
  if [[ $maxparallel -le 1 ]]; then
    checkDirectory $f
  else
    checkDirectory $f &
    job_limit $maxparallel
  fi
done

wait

echo "(OK:COPY_FAIL:FILE_DNE:FAIL:UNKNOWN) = (${nOK}:${nCOPYFAIL}:${nFILEDNE}:${nFAIL}:${nUNKNOWN})"
