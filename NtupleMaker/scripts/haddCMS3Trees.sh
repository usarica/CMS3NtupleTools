#!/bin/bash

JOBSDIR=""
OUTDIR=""
NCORES="12"
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
  exit 0
fi


for INDIR in $(ls ${JOBSDIR}); do
  for d in $(ls ${JOBSDIR}/${INDIR}); do
    pushd ${JOBSDIR}/${INDIR}/${d}
    RECDIR=${OUTDIR}/${INDIR}/${d}
    mkdir -p ${OUTDIR}/${INDIR}/${d}

    flist=""
    for f in $(ls ./ | grep -e ".root" | grep -v allevents.root); do
      if [[ -z "$flist" ]];then
        flist="$f"
      else
        flist="$flist $f"
      fi
    done

    hadd -ff -O -n 0 -j ${NCORES} ${RECDIR}/allevents.root $flist

    popd
  done
done
