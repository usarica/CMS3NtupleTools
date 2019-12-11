#!/bin/bash

INFILE=""
DATE=""
OUTPUTDIR=""
CONDOROUTDIR=""
QUEUE="vanilla"

let printhelp=0
for fargo in "$@";do
  fcnargname=""
  farg="${fargo//\"}"
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "infile="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    INFILE="$fcnargname"
  elif [[ "$fargl" == "date="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    DATE="$fcnargname"
  elif [[ "$fargl" == "outdir="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    OUTPUTDIR="$fcnargname"
  elif [[ "$fargl" == "condoroutdir="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    CONDOROUTDIR="$fcnargname"
  elif [[ "$fargl" == "help" ]];then
    let printhelp=1
  fi
done
if [[ $printhelp -eq 1 ]] || [[ -z "$INFILE" ]]; then
  echo "$0 usage:"
  echo " - help: Print this help"
  echo " - infile: Input commands list file. Mandatory."
  echo " - outdir: Main output location. Default='./output'"
  echo " - date: Date of the generation; does not have to be an actual date. Default=[today's date in YYMMDD format]"
  echo " - condoroutdir: Condor output directory to override. Optional."
  exit 0
fi

INITIALDIR=$(pwd)
if [[ "$FRAGDIR" != "/"* ]]; then
  FRAGDIR=${INITIALDIR}/${FRAGDIR}
fi
if [[ "$GRIDPACKDIR" != "/"* ]]; then
  GRIDPACKDIR=${INITIALDIR}/${GRIDPACKDIR}
fi

hname=$(hostname)

CONDORSITE="DUMMY"
if [[ "$hname" == *"lxplus"* ]];then
  echo "Setting default CONDORSITE to cern.ch"
  CONDORSITE="cern.ch"
elif [[ "$hname" == *"ucsd"* ]];then
  echo "Setting default CONDORSITE to t2.ucsd.edu"
  CONDORSITE="t2.ucsd.edu"
fi

if [[ "$OUTPUTDIR" == "" ]];then
  OUTPUTDIR="./output"
fi
if [[ "$DATE" == "" ]];then
  DATE=$(date +%y%m%d)
fi

OUTDIR="${OUTPUTDIR}/${DATE}"

mkdir -p $OUTDIR


TARFILE="cms3ntuplemaker.tar"
if [[ ! -e ${OUTDIR}/${TARFILE} ]];then
  createCMS3NtupleMakerTarball.sh
  mv ${TARFILE} ${OUTDIR}/
fi

checkGridProxy.sh

while IFS='' read -r line || [[ -n "$line" ]]; do
  THECONDORSITE="${CONDORSITE}"
  THECONDOROUTDIR="${CONDOROUTDIR}"
  THEQUEUE="${QUEUE}"
  THEOUTPUTFILE=""
  FCNARGS=""
  fcnarglist=($(echo $line))
  fcnargname=""
  for fargo in "${fcnarglist[@]}";do
    farg="${fargo//\"}"
    fargl="$(echo $farg | awk '{print tolower($0)}')"
    if [[ "$farg" == "output="* ]];then
      fcnargname=$farg
      fcnargname=${fcnargname//"output="}
      fcnargname=${fcnargname//".root"}
      THEOUTPUTFILE="${fcnargname}"
      echo "- Setting the output file to ${THEOUTPUTFILE}"
      if [[ -z "$FCNARGS" ]];then
        FCNARGS="$fargo"
      else
        FCNARGS="${FCNARGS} ${fargo}"
      fi
    elif [[ "$fargl" == "condorsite="* ]];then
      fcnargname="$farg"
      fcnargname="${fcnargname#*=}"
      THECONDORSITE="$fcnargname"
    elif [[ "$fargl" == "condoroutdir="* ]];then
      fcnargname="$farg"
      fcnargname="${fcnargname#*=}"
      THECONDOROUTDIR="$fcnargname"
    elif [[ "$fargl" == "condorqueue="* ]];then
      fcnargname="$farg"
      fcnargname="${fcnargname#*=}"
      THEQUEUE="$fcnargname"
    else
      if [[ -z "$FCNARGS" ]];then
        FCNARGS="$fargo"
      else
        FCNARGS="${FCNARGS} ${fargo}"
      fi
    fi
  done

  THECONDOROUTDIR=${THECONDOROUTDIR/"<USER>"/"$USER"}
  THECONDOROUTDIR=${THECONDOROUTDIR/"<DATE>"/"$DATE"}
  if [[ "${THECONDORSITE+x}" != "DUMMY" ]] && [[ -z "${THECONDOROUTDIR+x}" ]]; then
    echo "Need to set the Condor output directory."
    continue
  else
    echo "Condor directory chosen: ${THECONDORSITE}:${THECONDOROUTDIR}"
  fi

  theOutdir="${OUTDIR}/${THEOUTPUTFILE}"
  mkdir -p "${theOutdir}/Logs"

  ln -sf ${PWD}/${OUTDIR}/${TARFILE} ${PWD}/${theOutdir}/


  configureCMS3NtupleMakerCondorJob.py \
    --tarfile="$TARFILE" --batchqueue="$THEQUEUE" --outdir="$theOutdir" \
    --fcnargs="$FCNARGS" --condorsite="$THECONDORSITE" --condoroutdir="$THECONDOROUTDIR" \
    --outlog="Logs/log_job" --errlog="Logs/err_job" --batchscript="runCMS3NtupleMaker.condor.sh" --dry

done < "$INFILE"
