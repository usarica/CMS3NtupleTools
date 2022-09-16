#!/bin/bash

getvarpaths(){
  for var in "$@";do
    tmppath=${var//:/ }
    for p in $(echo $tmppath);do
      if [[ -e $p ]];then
        echo $p
      fi
    done
  done  
}
searchfileinvar(){
  for d in $(getvarpaths $1);do
    for f in $(ls $d | grep $2);do
      echo "$d/$f"
    done
  done
}
getcmssw(){
  if [ -r "$OSGVO_CMSSW_Path"/cmsset_default.sh ]; then
    echo "sourcing environment: source $OSGVO_CMSSW_Path/cmsset_default.sh"
    source "$OSGVO_CMSSW_Path"/cmsset_default.sh
  elif [ -r "$OSG_APP"/cmssoft/cms/cmsset_default.sh ]; then
    echo "sourcing environment: source $OSG_APP/cmssoft/cms/cmsset_default.sh"
    source "$OSG_APP"/cmssoft/cms/cmsset_default.sh
  elif [ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]; then
    echo "sourcing environment: source /cvmfs/cms.cern.ch/cmsset_default.sh"
    source /cvmfs/cms.cern.ch/cmsset_default.sh
  else
    echo "ERROR! Couldn't find $OSGVO_CMSSW_Path/cmsset_default.sh or /cvmfs/cms.cern.ch/cmsset_default.sh or $OSG_APP/cmssoft/cms/cmsset_default.sh"
    exit 1
  fi
}

CMSSWVERSION="$1"
SCRAMARCH="$2"
SUBMIT_DIR="$3" # Must be within $CMSSW_BASE/src/
TARFILE="$4"
SCRIPT="pebbleRoutine.cc"
FCN="pebbleRoutine"
FCNARGS=""

export SCRAM_ARCH=${SCRAMARCH}

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "SUBMIT_DIR: $SUBMIT_DIR"
echo "TARFILE: $TARFILE"
echo "SCRIPT: $SCRIPT"
echo "FCN: $FCN"
echo "FCNARGS: $FCNARGS"

echo "GLIDEIN_CMSSite: $GLIDEIN_CMSSite"
echo "hostname: $(hostname)"
echo "uname -a: $(uname -a)"
echo "whoami: $(whoami)"
echo "time: $(date +%s)"
echo "args: $@"
echo -e "\n--- end header output ---\n" #                       <----- section division

echo -e "\n--- begin memory specifications ---\n" #                     <----- section division
ulimit -a
echo -e "\n--- end memory specifications ---\n" #                     <----- section division


INITIALDIR=$(pwd)
echo "Initial directory is ${INITIALDIR}"

mkdir -p rundir
cd rundir

getcmssw

# If the first file in the tarball filelist starts with CMSSW, it is a
# tarball made outside of the full CMSSW directory and must be handled
# differently
if [[ ! -s ${INITIALDIR}/${TARFILE} ]];then
  echo "Tar file ${INITIALDIR}/${TARFILE} does not exist"
fi
if [[ ! -z $(tar -tf ${INITIALDIR}/${TARFILE} | head -n 1 | grep "^CMSSW") ]]; then

  echo "This is a full cmssw tar file."

  mv ${INITIALDIR}/${TARFILE} $(pwd)/
  if [[ "${TARFILE}" == *".tgz" ]];then
    tar zxf ${TARFILE}
  else
    tar xf ${TARFILE}
  fi
  rm ${TARFILE}

  if [[ -e extras.more ]];then
    mv extras.more ${CMSSWVERSION}/extras.tar
  fi

  cd $CMSSWVERSION

  if [[ -e extras.more ]];then
    tar xf extras.tar
    rm extras.tar
  fi

  echo "Current directory ${PWD} =? ${CMSSWVERSION}"
  echo "Running ProjectRename"
  scramv1 b ProjectRename

else

  # Setup the CMSSW area
  echo "This is a selective CMSSW tar file."
  eval `scramv1 project CMSSW $CMSSWVERSION`
  cd $CMSSWVERSION

  mv ${INITIALDIR}/${TARFILE} $(pwd)/
  if [[ "${TARFILE}" == *".tgz" ]];then
    tar zxf ${TARFILE}
  else
    tar xf ${TARFILE}
  fi
  rm ${TARFILE}

  if [[ -e extras.more ]];then
    mv extras.more extras.tar
    tar xf extras.tar
    rm extras.tar
  fi

fi


# Setup the CMSSW environment
eval `scramv1 runtime -sh`
echo "CMSSW_BASE: ${CMSSW_BASE}"

cd src

declare -i doRecompile
let doRecompile=0
if [[ -e CMS3/AnalysisTree/src ]];then
  declare -i nSrcFiles
  let nSrcFiles=$(ls CMS3/AnalysisTree/src | wc -l)
  if [[ $nSrcFiles -gt 0 ]];then
    let doRecompile=1
  fi
fi


# Set up environment variables for MELA and related packages
pushd ${CMSSW_BASE}/src &> /dev/null
eval $(./MelaAnalytics/setup.sh env)
popd &> /dev/null
# Needed to locate the include directory of MELA classes. It can get lost.
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${MELA_LIB_PATH}/../../interface


# Do the final compilation
if [[ $doRecompile -eq 1 ]]; then
  scramv1 b -j &>> compilation.log
else
  scramv1 b -j 1 &>> compilation.log
fi
CMSSW_COMPILE_STATUS=$?
if [ $CMSSW_COMPILE_STATUS != 0 ];then
  echo "CMSSW compilation exited with error ${CMSSW_COMPILE_STATUS}. Printing the log:"
  cat compilation.log
fi
rm -f compilation.log

# Go into the submission directory within $CMSSW_BASE/src
cd $SUBMIT_DIR

echo "Submission directory before running: ls -lrth"
ls -lrth

##############
# ACTUAL RUN #
##############
echo -e "\n--- Begin RUN ---\n"

# Create this file in order to detect in the scripts that a Condor job is running
touch RUNNING_ON_CONDOR

RUNDIR=$(pwd)
# Copy MELA-linked objects in case symlink does not work
## Taken from Mela.cc
cp -f ${CMSSW_BASE}/src/JHUGenMELA/MELA/data/process.DAT ./
cp -f ${CMSSW_BASE}/src/JHUGenMELA/MELA/data/br.sm1 ./
cp -f ${CMSSW_BASE}/src/JHUGenMELA/MELA/data/br.sm2 ./
cp -f ${CMSSW_BASE}/src/JHUGenMELA/MELA/data/ffwarn.dat ./
cp -f ${CMSSW_BASE}/src/JHUGenMELA/MELA/data/input.DAT ./
cp -rf ${CMSSW_BASE}/src/JHUGenMELA/MELA/data/Pdfdata ./
# Run the script
LOADLIB="loadLib.C"
RUN_CMD=$(runGenericROOTCommand.py --loadlib="$LOADLIB" --script="$SCRIPT" --function="$FCN" --command="$FCNARGS" --recompile --dry) # Must force recompilation
if [[ "$RUN_CMD" == "Running "* ]];then
  echo "$RUN_CMD"
  RUN_CMD=${RUN_CMD//"Running "/"timeout -s SIGKILL 10m "}
  eval "$RUN_CMD"
  RUN_STATUS=$?
  if [ $RUN_STATUS != 0 ]; then
    echo "Run has crashed with exit code ${RUN_STATUS}"
    sleep 1d
    exit 1
  else
    echo "Run succeeded. Pending for 20 min. to allow other pebbles to fill gaps..."
    sleep 20m
  fi
else
  echo "Run command ${RUN_CMD} is invalid."
  exit 1
fi

echo -e "\n--- End RUN ---\n"

echo "Submission directory after running all steps: ls -lrth"
ls -lrth

echo "time at end: $(date +%s)"
