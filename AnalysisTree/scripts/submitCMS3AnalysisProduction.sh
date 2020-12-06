#!/bin/bash

SCRIPT=""
FCN=""
FCNARGS=""
DATE=""
OUTPUTDIR=""
CONDOROUTDIR="/hadoop/cms/store/user/<USER>/Offshell_2L2Nu/Worker"
EXTRATARCMD=""
QUEUE="vanilla"
REQMEM="2048M"
JOBFLAVOR="tomorrow"

let recreate=0
let checkproxy=1
let printhelp=0
for fargo in "$@";do
  fcnargname=""
  farg="${fargo//\"}"
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "script="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    SCRIPT="$fcnargname"
  elif [[ "$fargl" == "function="* ]] || [[ "$fargl" == "fcn="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    FCN="$fcnargname"
  elif [[ "$fargl" == "arguments="* ]] || [[ "$fargl" == "fcnargs="* ]];then
    fcnargname="$fargo"
    fcnargname="${fcnargname#*=}"
    FCNARGS="$fcnargname"
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
  elif [[ "$fargl" == "memory="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    REQMEM="$fcnargname"
  elif [[ "$fargl" == "job_flavor="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    JOBFLAVOR="$fcnargname"
  elif [[ "$fargl" == "tarinclude="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    if [[ -z "$EXTRATARCMD" ]];then
      EXTRATARCMD="addfile=$fcnargname"
    else
      EXTRATARCMD="$EXTRATARCMD addfile=$fcnargname"
    fi
  elif [[ "$fargl" == "recreate" ]];then
    let recreate=1
  elif [[ "$fargl" == "no-proxycheck" ]];then
    let checkproxy=0
  elif [[ "$fargl" == "help" ]];then
    let printhelp=1
  fi
done
if [[ $printhelp -eq 1 ]] || [[ -z "$SCRIPT" ]]; then
  echo "$0 usage:"
  echo " - help: Print this help"
  echo " - script: Name of the script. Mandatory."
  echo " - function / fcn: Name of the function in the script. Default=[script]"
  echo " - arguments / fcnargs: Arguments of the function. Default=\"\""
  echo " - outdir: Main output location. Default='./output'"
  echo " - date: Date of the generation; does not have to be an actual date. Default=[today's date in YYMMDD format]"
  echo " - condoroutdir: Condor output directory to override. Default=/hadoop/cms/store/user/<USER>/Offshell_2L2Nu/Worker"
  echo " - memory: Required RAM for the job. Default='2048M'"
  echo " - job_flavor: Required time limit flavor for the job. Can be 'workday', 'tomorrow', 'testmatch', or 'nextweek'. Default='tomorrow'"
  echo " - tarinclude: Include extra files in the tar. Can specify multiple times. Default: None"
  echo " - recreate: Force the recreation of the job directories"
  echo " - no-proxycheck: Do not check the proxy"
  exit 0
fi
SCRIPTRAWNAME="${SCRIPT%%.*}"
if [[ -z "$FCN" ]];then
  FCN=$SCRIPTRAWNAME
fi

INITIALDIR=$(pwd)

hname=$(hostname)

if [[ -z "$OUTPUTDIR" ]];then
  OUTPUTDIR="./output"
fi
if [[ -z "$DATE" ]];then
  DATE=$(date +%y%m%d)
fi

OUTDIR="${OUTPUTDIR}/${DATE}"

mkdir -p $OUTDIR

TARFILE="cms3analysistree.tar"
if [[ ! -e ${OUTDIR}/${TARFILE} ]];then
  createCMS3AnalysisTreeTarball.sh ${EXTRATARCMD}
  mv ${TARFILE} ${OUTDIR}/
fi

THEOUTPUTFILE=$SCRIPTRAWNAME
if [[ "$FCN" != "$SCRIPTRAWNAME" ]];then
  THEOUTPUTFILE="${THEOUTPUTFILE}/${FCN}"
fi
if [[ ! -z "${FCNARGS}" ]];then
  fcnargname="${FCNARGS}"
  fcnargname="${fcnargname// /}"
  fcnargname="${fcnargname//,/_}"
  fcnargname="${fcnargname//=/_}"
  fcnargname="${fcnargname//.root}"
  fcnargname="${fcnargname//\"}"
  fcnargname="${fcnargname//\!}"
  fcnargname="${fcnargname//\\}"
  fcnargname="${fcnargname//(}"
  fcnargname="${fcnargname//)}"
  fcnargname="${fcnargname//./p}"
  if [[ "$fcnargname" == "/"* ]];then
    fcnargname="${fcnargname:1}"
  fi
  fcnargname="${fcnargname//\//_}"
  THEOUTPUTFILE="${THEOUTPUTFILE}/${fcnargname}"
else
  FCNARGS="<NONE>"
fi

theOutdir="${OUTDIR}/${THEOUTPUTFILE}"
let hasJobSetup=0
if [[ -e $theOutdir ]] || [[ -e ${theOutdir}.tar ]]; then
  let hasJobSetup=1
fi
if [[ $recreate -eq 1 ]] || [[ $hasJobSetup -eq 0 ]]; then
  CONDORSITE="DUMMY"
  if [[ "$hname" == *"lxplus"* ]];then
    echo "Setting default CONDORSITE to cern.ch"
    CONDORSITE="cern.ch"
  elif [[ "$hname" == *"ucsd"* ]];then
    echo "Setting default CONDORSITE to t2.ucsd.edu"
    CONDORSITE="t2.ucsd.edu"
  fi

  THECONDORSITE="${CONDORSITE}"
  THECONDOROUTDIR="${CONDOROUTDIR}"
  THEQUEUE="${QUEUE}"
  THECONDOROUTDIR=${THECONDOROUTDIR/"<USER>"/"$USER"}
  THECONDOROUTDIR=${THECONDOROUTDIR/"<DATE>"/"$DATE"}
  THECONDOROUTDIR=${THECONDOROUTDIR/"<SCRIPT>"/"$SCRIPTRAWNAME"}
  THECONDOROUTDIR=${THECONDOROUTDIR/"<FCN>"/"$FCN"}
  if [[ "${THECONDORSITE+x}" != "DUMMY" ]] && [[ -z "${THECONDOROUTDIR+x}" ]]; then
    echo "Need to set the Condor output directory."
    continue
  else
    echo "Condor directory chosen: ${THECONDORSITE}:${THECONDOROUTDIR}"
  fi

  if [[ $checkproxy -eq 1 ]]; then
    checkGridProxy.sh
  fi

  mkdir -p "${theOutdir}/Logs"

  ln -sf ${PWD}/${OUTDIR}/${TARFILE} ${PWD}/${theOutdir}/

  configureCMS3AnalysisTreeCondorJob.py \
    --tarfile="$TARFILE" --batchqueue="$THEQUEUE" --outdir="$theOutdir" \
    --script="$SCRIPT" --fcn="$FCN" --fcnargs="${FCNARGS}" --required_memory="${REQMEM}" --job_flavor="${JOBFLAVOR}" \
    --condorsite="$THECONDORSITE" --condoroutdir="$THECONDOROUTDIR" \
    --outlog="Logs/log_job" --errlog="Logs/err_job" --batchscript="runCMS3AnalysisTree.condor.sh" --dry
fi
