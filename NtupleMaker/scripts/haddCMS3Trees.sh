#!/bin/bash

JOBSDIR=""
OUTDIR=""
NCORES="12"
let OVERWRITE=0
let OPTIMIZE=0
let printhelp=0
for fargo in "$@";do
  fcnargname=""
  farg="${fargo//\"}"
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "jobsdir="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    JOBSDIR="$fcnargname"
  elif [[ "$fargl" == "outdir="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    OUTDIR="$fcnargname"
  elif [[ "$fargl" == "ncores="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    NCORES="$fcnargname"
  elif [[ "$fargl" == "overwrite="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    fcnargname="$(echo $fcnargname | awk '{print tolower($0)}')"
    if [[ -z "$fcnargname" ]] || [[ "$fcnargname" == "t"* ]] || [[ "$fcnargname" == "y"* ]] || [[ "$fcnargname" == "1" ]];then
      let OVERWRITE=1
    elif [[ "$fcnargname" == "f"* ]] || [[ "$fcnargname" == "n"* ]] || [[ "$fcnargname" == "0" ]];then
      let OVERWRITE=0
    else
      let printhelp=1
    fi
  elif [[ "$fargl" == "optimize="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    fcnargname="$(echo $fcnargname | awk '{print tolower($0)}')"
    if [[ -z "$fcnargname" ]] || [[ "$fcnargname" == "t"* ]] || [[ "$fcnargname" == "y"* ]] || [[ "$fcnargname" == "1" ]];then
      let OPTIMIZE=1
    elif [[ "$fcnargname" == "f"* ]] || [[ "$fcnargname" == "n"* ]] || [[ "$fcnargname" == "0" ]];then
      let OPTIMIZE=0
    else
      let printhelp=1
    fi
  elif [[ "$fargl" == "help" ]];then
    let printhelp=1
  fi
done
if [[ $printhelp -eq 1 ]] || [[ -z "$JOBSDIR" ]] || [[ -z "$OUTDIR" ]]; then
  echo "$0 usage:"
  echo " - help: Print this help"
  echo " - jobsdir: Job submission directory. Mandatory."
  echo " - outdir: Main output location. Mandatory."
  echo " - ncores: Number of cores for hadd (Default: 12)"
  echo " - overwrite: Flag to specify whether the output file should be overwritten. For true value, specify either nothing after the option, or t*/y*/1. For false value, either do not specify this option, or put f*/n*/0. (Default: false)"
  echo " - optimize: Flag to optimize file size. (Default: false)"
  exit 0
fi


for INDIR in $(ls ${JOBSDIR}); do
  for d in $(ls ${JOBSDIR}/${INDIR}); do
    JOBSINDIR=${JOBSDIR}/${INDIR}/${d}
    RECDIR=${OUTDIR}/${INDIR}/${d}
    OUTFILE=${RECDIR}/allevents.root

    if [[ $OVERWRITE -eq 0 ]] && [[ -s ${OUTFILE} ]];then
      echo "${OUTFILE} exists. Skipping..."
      continue
    fi

    mkdir -p ${RECDIR}

    pushd ${JOBSINDIR}

    flist=""
    for f in $(ls ./ | grep -e ".root" | grep -v allevents.root); do
      if [[ -z "$flist" ]];then
        flist="$f"
      else
        flist="$flist $f"
      fi
    done

    optOptimize=""
    if [[ ${OPTIMIZE} -eq 1 ]];then
      optOptimize="-O"
    fi

    hadd -ff ${optOptimize} -n 0 -j ${NCORES} ${OUTFILE} $flist

    popd
  done
done
