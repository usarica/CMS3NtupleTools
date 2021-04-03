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
TARFILE="$3"
WSTARFILE="$4"
CONDORSITE="$5"
CONDOROUTDIR="$6"


export SCRAM_ARCH=${SCRAMARCH}

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "TARFILE: $TARFILE"
echo "WSTARFILE: $WSTARFILE"
echo "CONDORSITE: $CONDORSITE"
echo "CONDOROUTDIR: $CONDOROUTDIR"

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

# Do the final compilation
scramv1 b -j 1 &>> compilation.log
CMSSW_COMPILE_STATUS=$?
if [ $CMSSW_COMPILE_STATUS != 0 ];then
  echo "CMSSW compilation exited with error ${CMSSW_COMPILE_STATUS}. Printing the log:"
  cat compilation.log
  exit ${CMSSW_COMPILE_STATUS}
fi
rm -f compilation.log

# Go into the submission directory within $CMSSW_BASE/src
mkdir -p subdir
cd subdir

mv ${INITIALDIR}/${WSTARFILE} ./
if [[ "${WSTARFILE}" == *".tgz" ]];then
  tar zxf ${WSTARFILE}
else
  tar xf ${WSTARFILE}
fi
rm ${WSTARFILE}

# The data card and workspace file might be in a subdirectory, so we need to move them here.
DCFILE=$(find ./ -name datacard.txt)
WSFILE=$(find ./ -name workspace.root)
mv ${DCFILE} ./
mv ${WSFILE} ./

echo "Submission directory before running: ls -lrth"
ls -lrth

##############
# ACTUAL RUN #
##############
echo -e "\n--- Begin RUN ---\n"

# Create this file in order to detect in the scripts that a Condor job is running
touch RUNNING_ON_CONDOR

RUNDIR=$(pwd)

# Assign the interactive library loading instructions from Combine
LOADLIB="${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/loadLib.C"

# Compile the list of commands
CMDSFILE=cmds.txt
echo "text2workspace.py --X-allow-no-background datacard.txt -o combined.root" >> ${CMDSFILE}
FITOPTS="-M MultiDimFit --freezeParameters r --redefineSignalPOIs frac_sig -v 3 --alignEdges=1 --saveNLL --saveSpecifiedNuis=rgx{.*} --algo singles" # --saveFitResult
echo "combine combined.root ${FITOPTS} -n SinglesSnapshot1 --saveWorkspace" >> ${CMDSFILE}
echo "combine higgsCombineSinglesSnapshot1.MultiDimFit.mH120.root ${FITOPTS} -n SinglesSnapshot2 --snapshotName=MultiDimFit --saveWorkspace" >> ${CMDSFILE}
echo "combine higgsCombineSinglesSnapshot2.MultiDimFit.mH120.root ${FITOPTS} -n SinglesFinalSnapshot --snapshotName=MultiDimFit --saveWorkspace" >> ${CMDSFILE}
echo "mv higgsCombineSinglesFinalSnapshot.MultiDimFit.mH120.root combined_withSnapshot.root" >> ${CMDSFILE}

# Run the commands
while IFS='' read -r RUN_CMD || [[ -n "$RUN_CMD" ]]; do
  echo "Running ${RUN_CMD}"
  eval "$RUN_CMD"
  RUN_STATUS=$?
  if [ $RUN_STATUS != 0 ]; then
    echo "Run has crashed with exit code ${RUN_STATUS}"
    exit 1
  fi
done < "${CMDSFILE}"

echo -e "\n--- End RUN ---\n"


echo "Submission directory after running all steps before file transfers: ls -lrth"
ls -lrth


##################
# TRANSFER FILES #
##################
echo -e "\n--- Begin transfer ---\n"

# Once you reach here, you have the file 'combined_withSnapshot.root'
OUTFILENAME="combined_withSnapshot.root"
echo "Copying output file ${OUTFILENAME}"
# Use the executable 'copyFromCondorToSite.sh' from CMSDataTools/AnalysisTree already in ${CMSSW_BASE}/bin/.
copyFromCondorToSite.sh ${RUNDIR} ${OUTFILENAME} ${CONDORSITE} ${CONDOROUTDIR}
TRANSFER_STATUS=$?
if [ $TRANSFER_STATUS != 0 ]; then
  echo " - Transfer crashed with exit code ${TRANSFER_STATUS}"
  exit ${TRANSFER_STATUS}
fi

echo -e "\n--- End transfer ---\n"


echo "Submission directory after running all steps: ls -lrth"
ls -lrth

echo "time at end: $(date +%s)"
