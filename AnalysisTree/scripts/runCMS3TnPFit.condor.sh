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
declare -i RERUNWITHROBUSTFIT=0
for fargo in "$@"; do
  farg="${fargo//\"}"
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "rerun_with_robustfit" ]]; then
    RERUNWITHROBUSTFIT=1
  fi
done

export SCRAM_ARCH=${SCRAMARCH}

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "TARFILE: $TARFILE"
echo "WSTARFILE: $WSTARFILE"
echo "CONDORSITE: $CONDORSITE"
echo "CONDOROUTDIR: $CONDOROUTDIR"
echo "RERUNWITHROBUSTFIT: $RERUNWITHROBUSTFIT"

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

# Make sure JHUGenMELA library location is recognized
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CMSSW_BASE}/src/JHUGenMELA/MELA/data/${SCRAM_ARCH}

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

if [[ -f ${INITIALDIR}/${WSTARFILE} ]]; then
  mv ${INITIALDIR}/${WSTARFILE} ./
  if [[ "${WSTARFILE}" == *".tgz" ]];then
    tar zxf ${WSTARFILE}
  else
    tar xf ${WSTARFILE}
  fi
  rm ${WSTARFILE}
fi

# The data card and workspace file might be in a subdirectory, so we need to move them here.
DCFILE=$(find ./ -name datacard.txt)
WSFILE=$(find ./ -name workspace.root)
if [[ "${DCFILE}" != "" ]]; then
  mv ${DCFILE} ./
elif [[ ${RERUNWITHROBUSTFIT} -eq 0 ]]; then
  echo "Cannot find datacard.txt"
  exit 1
fi
if [[ "${WSFILE}" != "" ]]; then
  mv ${WSFILE} ./
elif [[ ${RERUNWITHROBUSTFIT} -eq 0 ]]; then
  echo "Cannot find workspace.root"
  exit 1
fi

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
FITOPTS=""
declare -a output_files=( )
if [[ ${RERUNWITHROBUSTFIT} -eq 0 ]]; then
  # Configure the main fit
  FITOPTS="-M MultiDimFit --freezeParameters r --redefineSignalPOIs frac_sig -v 3 --alignEdges=1 --saveNLL --saveSpecifiedNuis=rgx{.*} --algo singles --saveWorkspace --X-rtd MINIMIZER_no_analytic"
  finalOutputFile=combined_withSnapshot.root

  echo "timeout -s SIGKILL 5m text2workspace.py --no-b-only --X-allow-no-background datacard.txt -o combined.root" >> ${CMDSFILE}
  echo "combine combined.root ${FITOPTS} -n SinglesSnapshot1" >> ${CMDSFILE}
  echo "combine higgsCombineSinglesSnapshot1.MultiDimFit.mH120.root ${FITOPTS} -n SinglesSnapshot2 --snapshotName=MultiDimFit" >> ${CMDSFILE}
  echo "combine higgsCombineSinglesSnapshot2.MultiDimFit.mH120.root ${FITOPTS} -n SinglesFinalSnapshot --snapshotName=MultiDimFit" >> ${CMDSFILE}
  echo "mv higgsCombineSinglesFinalSnapshot.MultiDimFit.mH120.root ${finalOutputFile}" >> ${CMDSFILE}
  output_files+=( $finalOutputFile )

  # Configure the 1-sigma scan using robust fit
  FITOPTS="${FITOPTS} --robustFit 1 --X-rtd FITTER_CROSSING_PERSISTENT_TRIALS=30"
  initialWSName=${finalOutputFile}
  finalOutputFile=combined_withSnapshot_withRobustFit.root

  echo "combine ${initialWSName} ${FITOPTS} -n SinglesFinalSnapshotWithRobustFit --snapshotName=MultiDimFit" >> ${CMDSFILE}
  echo "mv higgsCombineSinglesFinalSnapshotWithRobustFit.MultiDimFit.mH120.root ${finalOutputFile}" >> ${CMDSFILE}

  output_files+=( $finalOutputFile )
else
  FITOPTS="-M MultiDimFit --freezeParameters r --redefineSignalPOIs frac_sig -v 3 --alignEdges=1 --saveNLL --saveSpecifiedNuis=rgx{.*} --algo singles --saveWorkspace --X-rtd MINIMIZER_no_analytic --robustFit 1 --X-rtd FITTER_CROSSING_PERSISTENT_TRIALS=30"
  initialWSName="$(ExecuteCompiledCommand GetStandardHostPathToStore ${CONDOROUTDIR}/combined_withSnapshot.root ${CONDORSITE})"
  finalOutputFile=combined_withSnapshot_withRobustFit.root

  echo "combine ${initialWSName} ${FITOPTS} -n SinglesFinalSnapshotWithRobustFit --snapshotName=MultiDimFit" >> ${CMDSFILE}
  echo "mv higgsCombineSinglesFinalSnapshotWithRobustFit.MultiDimFit.mH120.root ${finalOutputFile}" >> ${CMDSFILE}

  output_files+=( $finalOutputFile )
fi

if [[ -f datacard.txt ]]; then
  strNData="$(grep rate datacard.txt | awk '{print $2}')"
  if [[ "$strNData" == "0" ]]; then
    echo "Datacard has data channel with 0 events:"
    grep rate datacard.txt
    echo "ZERO_DATA_EVENT" >> VALID_FAIL.txt
    output_files=( VALID_FAIL.txt )
    rm -f ${CMDSFILE}
    touch ${CMDSFILE}
  fi
  strRateLine="$(grep rate datacard.txt)"
  strRateLineNew="${strRateLine//0.00000/0.00001}"
  sed -i "s|${strRateLine}|${strRateLineNew}|g" datacard.txt

  echo "Data card to fit:"
  cat datacard.txt
fi

echo "Will execute the following commands:"
cat ${CMDSFILE}

# Run the commands
while IFS='' read -r RUN_CMD || [[ -n "$RUN_CMD" ]]; do
  rm -f runlog.log
  rm -f rerunlog.log
  echo "Running ${RUN_CMD}"
  eval "${RUN_CMD}" |& tee runlog.log
  RUN_STATUS=$?

  RERUN_TWS=0
  RERUN_COMBINE=0
  if [[ ${RUN_STATUS} -ne 0 ]] && [[ "${RUN_CMD}" == *"text2workspace"* ]]; then
    let RERUN_TWS=1
  elif [[ "${RUN_CMD}" == "combine"* ]]; then
    failstr="$(grep -e 'WARNING: MultiDimFit failed' runlog.log)"
    if [[ ${RUN_STATUS} -ne 0 ]] || [[ "${failstr}" != "" ]]; then
      let RERUN_COMBINE=1
    fi
  fi

  if [[ ${RERUN_TWS} -eq 1 ]]; then
    rm -f combined.root
    trimCMS3TnPWSDatasets workspace.root size_req=62914560
    RERUN_CMD="${RUN_CMD/datacard.txt/datacard_trimmedTnPData.txt}"
    echo "Running ${RERUN_CMD} after trimming the workspace data set..."
    eval "$RERUN_CMD"
    RUN_STATUS=$?
  fi

  if [[ ${RERUN_COMBINE} -eq 1 ]]; then
    failedInput="$(echo ${RUN_CMD} | awk '{print $2}')"
    failedOutput="$(echo ${RUN_CMD##*-n })"
    failedOutput="$(echo ${failedOutput%% *})"
    echo "Recovering failed fit in ${failedOutput}"
    RERUN_CMD="${RUN_CMD} --cminDefaultMinimizerStrategy 0"
    RERUN_CMD="${RERUN_CMD/${failedInput}/higgsCombine${failedOutput}.MultiDimFit.mH120.root}"
    RERUN_CMD="${RERUN_CMD/-n ${failedOutput}/-n ${failedOutput}_recovery}"
    echo "Running ${RERUN_CMD}"
    eval "$RERUN_CMD"
    RUN_STATUS=$?
    if [[ ${RUN_STATUS} -eq 0 ]]; then
      mv higgsCombine${failedOutput}_recovery.MultiDimFit.mH120.root higgsCombine${failedOutput}.MultiDimFit.mH120.root
    else
      echo "Fit recovery run failed with exit code ${RUN_STATUS}"
      rm higgsCombine${failedOutput}_recovery.MultiDimFit.mH120.root
      let RUN_STATUS=0
    fi
  fi

  if [[ ${RUN_STATUS} -ne 0 ]]; then
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

for OUTFILENAME in "${output_files[@]}"; do
  echo "Copying output file ${OUTFILENAME}"
  # Use the executable 'copyFromCondorToSite.sh' from CMSDataTools/AnalysisTree already in ${CMSSW_BASE}/bin/.
  copyFromCondorToSite.sh ${RUNDIR} ${OUTFILENAME} ${CONDORSITE} ${CONDOROUTDIR//.tar}
  TRANSFER_STATUS=$?
  if [ $TRANSFER_STATUS != 0 ]; then
    echo " - Transfer crashed with exit code ${TRANSFER_STATUS}"
    exit ${TRANSFER_STATUS}
  fi
done

echo -e "\n--- End transfer ---\n"


echo "Submission directory after running all steps: ls -lrth"
ls -lrth

echo "time at end: $(date +%s)"
