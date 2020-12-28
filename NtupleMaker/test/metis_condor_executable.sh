#!/bin/bash

OUTPUTDIR=$1
OUTPUTNAME=$2
INPUTFILENAMES=$3
IFILE=$4
PSET=$5
CMSSWVERSION=$6
SCRAMARCH=$7
NEVTS=$8
FIRSTEVT=$9
EXPECTEDNEVTS=${10}
OTHEROUTPUTS=${11}
PSETARGS="${@:12}" # since args can have spaces, we take 10th-->last argument as one

# Make sure OUTPUTNAME doesn't have .root since we add it manually
OUTPUTNAME=$(echo $OUTPUTNAME | sed 's/\.root//')

export SCRAM_ARCH=${SCRAMARCH}

function getjobad {
    grep -i "^$1" "$_CONDOR_JOB_AD" | cut -d= -f2- | xargs echo
}
function setup_chirp {
    if [ -e ./condor_chirp ]; then
        # Note, in the home directory
        mkdir chirpdir
        mv condor_chirp chirpdir/
        export PATH="${PATH}:$(pwd)/chirpdir"
        echo "[chirp] Found and put condor_chirp into $(pwd)/chirpdir"
    elif [ -e /usr/libexec/condor/condor_chirp ]; then
        export PATH="${PATH}:/usr/libexec/condor"
        echo "[chirp] Found condor_chirp in /usr/libexec/condor"
    else
        echo "[chirp] No condor_chirp :("
    fi
}
function chirp {
    command -v condor_chirp &> /dev/null
    if [[ $? -ne 0 ]]; then
      exit 0 # Just exit normally
    fi
    # Note, $1 (the classad name) must start with Chirp
    condor_chirp set_job_attr_delayed $1 $2
    ret=$?
    echo "[chirp] Chirped $1 => $2 with exit code $ret"
}

setup_chirp

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "OUTPUTDIR: $OUTPUTDIR"
echo "OUTPUTNAME: $OUTPUTNAME"
echo "INPUTFILENAMES: $INPUTFILENAMES"
echo "IFILE: $IFILE"
echo "PSET: $PSET"
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "NEVTS: $NEVTS"
echo "EXPECTEDNEVTS: $EXPECTEDNEVTS"
echo "OTHEROUTPUTS: $OTHEROUTPUTS"
echo "PSETARGS: $PSETARGS"
# echo  CLASSAD: $(cat "$_CONDOR_JOB_AD")

if [[ ! -z ${GLIDEIN_CMSSite+x} ]]; then
  echo "GLIDEIN_CMSSite: $GLIDEIN_CMSSite"
fi
echo "hostname: $(hostname)"
echo "uname -a: $(uname -a)"
echo "time: $(date +%s)"
echo "args: $@"
echo "tag: $(getjobad tag)"
echo "taskname: $(getjobad taskname)"

echo -e "\n--- end header output ---\n" #                       <----- section division

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

# holy crap this is a mess. :( why does PAT code have to do such insane
# things with paths?
# if the first file in the tarball filelist starts with CMSSW, then it is
# a tarball made outside of the full CMSSW directory, and must be handled
# differently
tarfile=package.tar.gz
if [ ! -z $(tar -tf ${tarfile} | head -n 1 | grep "^CMSSW") ]; then
    echo "this is a full cmssw tar file"
    tar xf ${tarfile}
    cd $CMSSWVERSION
    echo $PWD
    echo "Running ProjectRename"
    scramv1 b ProjectRename
    echo "Running `scramv1 runtime -sh`"
    eval `scramv1 runtime -sh`
    mv ../$PSET pset.py
    mv ../${tarfile} .
else
    echo "this is a selective cmssw tar file"
    eval `scramv1 project CMSSW $CMSSWVERSION`
    cd $CMSSWVERSION
    eval `scramv1 runtime -sh`
    mv ../$PSET pset.py
    if [ -e ../${tarfile} ]; then
        mv ../${tarfile} ${tarfile}
        tar xf ${tarfile}
    fi
    scram b
fi

# Needed or else cmssw can't find libmcfm_[xyz].so
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CMSSW_BASE}/src/JHUGenMELA/MELA/data/${SCRAM_ARCH}
# This is nicer than above. both work, and both have scary but benign warnings/printouts
cp ${CMSSW_BASE}/src/JHUGenMELA/MELA/data/${SCRAM_ARCH}/*.so ${CMSSW_BASE}/lib/${SCRAM_ARCH}/
# "Needed" to get rid of benign warnings/printouts
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${CMSSW_BASE}/src/JHUGenMELA/MELA/interface


# # logging every 45 seconds gives ~100kb log file/3 hours
# dstat -cdngytlmrs --float --nocolor -T --output dsout.csv 180 >& /dev/null &


echo "before running: ls -lrth"
ls -lrth

echo -e "\n--- begin running ---\n" #                           <----- section division

chirp ChirpMetisExpectedNevents $EXPECTEDNEVTS

chirp ChirpMetisStatus "before_cmsRun"

cmdRun="cmsRun pset.py inputs=${INPUTFILENAMES} output=${OUTPUTNAME}.root ${PSETARGS}"
echo "Running: ${cmdRun}"
${cmdRun}
CMSRUN_STATUS=$?

chirp ChirpMetisStatus "after_cmsRun"

echo "after running: ls -lrth"
ls -lrth

if [[ $CMSRUN_STATUS != 0 ]]; then
    echo "Removing output file because cmsRun crashed with exit code $?"
    rm ${OUTPUTNAME}.root
    exit 1
fi

CheckFileIntegrity ${OUTPUTNAME}.root
INTEGRITY_STATUS=$?
if [[ $INTEGRITY_STATUS != 0 ]]; then
    echo "Removing output file because integrity check failed with exit code $?"
    rm ${OUTPUTNAME}.root
    exit 1
fi

echo -e "\n--- end running ---\n" #                             <----- section division

echo "Sending output file ${OUTPUTNAME}.root"

if [ ! -e "${OUTPUTNAME}.root" ]; then
    echo "ERROR! Output ${OUTPUTNAME}.root doesn't exist"
    exit 1
fi

echo "time before copy: $(date +%s)"
chirp ChirpMetisStatus "before_copy"

copyFromCondorToSite.sh $(pwd) ${OUTPUTNAME}.root t2.ucsd.edu ${OUTPUTDIR} ${OUTPUTNAME}_${IFILE}.root

for OTHEROUTPUT in $(echo "$OTHEROUTPUTS" | sed -n 1'p' | tr ',' '\n'); do
    [ -e ${OTHEROUTPUT} ] && {
        NOROOT=$(echo $OTHEROUTPUT | sed 's/\.root//')
        copyFromCondorToSite.sh $(pwd) ${NOROOT}.root t2.ucsd.edu ${OUTPUTDIR} ${NOROOT}_${IFILE}.root
    }
done

echo -e "\n--- begin dstat output ---\n" #                      <----- section division
# cat dsout.csv
echo -e "\n--- end dstat output ---\n" #                        <----- section division

echo "time at end: $(date +%s)"

chirp ChirpMetisStatus "done"

