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
FCNARGS="$5"
CONDORSITE="$6"
CONDOROUTDIR="$7"

export SCRAM_ARCH=${SCRAMARCH}

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "SUBMIT_DIR: $SUBMIT_DIR"
echo "TARFILE: $TARFILE"
echo "FCNARGS: $FCNARGS"
echo "CONDORSITE: $CONDORSITE"
echo "CONDOROUTDIR: $CONDOROUTDIR"

echo "GLIDEIN_CMSSite: $GLIDEIN_CMSSite"
echo "hostname: $(hostname)"
echo "uname -a: $(uname -a)"
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
  tar xf ${TARFILE}
  rm ${TARFILE}

  cd $CMSSWVERSION
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
  tar xf ${TARFILE}
  rm ${TARFILE}
fi


# Setup the CMSSW environment
eval `scramv1 runtime -sh`
echo "CMSSW_BASE: ${CMSSW_BASE}"

cd src

# Clean CMSSW-related compilation objects and print the lib area afterward
scramv1 b clean &>> compilation.log
echo "================================="
echo "lib/${SCRAM_ARCH} after cleaning:"
ls ../lib/${SCRAM_ARCH}
echo "================================="


# Needed to locate the include directory of MELA classes. It can get lost.
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${CMSSW_BASE}/src/ZZMatrixElement/MELA/interface
# Ensure CMSSW can find libmcfm
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/${SCRAM_ARCH}
# Do not do the one below instead of the above; it will create problems when loading the MELA library interactively
# cp ${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/${SCRAM_ARCH}/*.so ${CMSSW_BASE}/lib/${SCRAM_ARCH}/


# Compile CMSSW-dependent packages
(
  cd ZZMatrixElement

  ./setup.sh clean
  ./setup.sh -j &>> compilation.log

  MELA_COMPILE_STATUS=$?
  if [ $MELA_COMPILE_STATUS != 0 ];then
    echo "MELA compilation exited with error ${MELA_COMPILE_STATUS}. Printing the log:"
    cat compilation.log
  fi
  rm -f compilation.log

  cd -
)
scramv1 b -j &>> compilation.log
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

RUNDIR=$(pwd)
# Copy MELA-linked objects in case symlink does not work
## Taken from Mela.cc
cp -f ${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/process.DAT ./
cp -f ${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/br.sm1 ./
cp -f ${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/br.sm2 ./
cp -f ${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/ffwarn.dat ./
cp -f ${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/input.DAT ./
cp -rf ${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/Pdfdata ./
# Run main_pset.py
RUNFILE=main_pset.py
RUN_CMD=$(runGenericExecutable.py --executable="$RUNFILE" --command="$FCNARGS" --dry)
if [[ "$RUN_CMD" == "Running "* ]];then
  echo "$RUN_CMD"
  RUN_CMD=${RUN_CMD//"Running "}
  eval "cmsRun -n 1 $RUN_CMD"
  RUN_STATUS=$?
  if [ $RUN_STATUS != 0 ]; then
    echo "Run has crashed with exit code ${RUN_STATUS}"
    exit 1
  fi
else
  echo "Run command ${RUN_CMD} is invalid."
  exit 1
fi

echo -e "\n--- End RUN ---\n"


##################
# TRANSFER FILES #
##################
echo "Submission directory after running all steps: ls -lrth"
ls -lrth

OUTFILENAME=""
fcnarglist=($(echo $RUN_CMD))
for fargo in "${fcnarglist[@]}";do
  farg="${fargo//\"}"
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "output="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    #fcnargname="${fcnargname//.root}"
    fcnargname="${fcnargname##*/}"
    OUTFILENAME="${fcnargname}"
  fi
done

if [[ ! -z ${OUTFILENAME} ]];then
  echo "Copying output file ${OUTFILENAME}"
  copyFromCondorToSite.sh ${RUNDIR} ${OUTFILENAME} ${CONDORSITE} ${CONDOROUTDIR}
fi
##############


echo "Submission directory after running: ls -lrth"
ls -lrth

echo "time at end: $(date +%s)"
