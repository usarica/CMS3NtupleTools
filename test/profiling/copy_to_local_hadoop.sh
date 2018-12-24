#!/usr/bin/env bash

# Because xrootd is a piece of shit when reading files piece by piece,
# copy over the whole file to hadoop first, and then simply use a slightly
# modified path to speed up reading and decrease the probability of violence

function docache {
    # Copies a file like /store/data/Run2017F/DoubleEG/MINIAOD/09May2018-v1/10000/444E03EB-B75F-E811-AFBA-F01FAFD8F16A.root
    # to /hadoop/cms/store/user/${USER}/localcache/data/Run2017F/DoubleEG/MINIAOD/09May2018-v1/10000/444E03EB-B75F-E811-AFBA-F01FAFD8F16A.root
    # I.e., /store/ is replaced with /hadoop/..../localcache/
    dest=$(echo "$1" | sed 's|/store|/hadoop/cms/store/user/'"$USER"'/localcache|')
    destdir=$(dirname $dest)
    bname=$(basename "$1")
    echo "xrdcp -S 4 root://cmsxrootd.fnal.gov/$1 $bname && mkdir -p $destdir && mv $bname $dest"
    xrdcp -S 4 root://cmsxrootd.fnal.gov/$1 $bname && mkdir -p $destdir && mv $bname $dest
    echo $dest
}

# docache /store/data/Run2017F/DoubleEG/MINIAOD/09May2018-v1/10000/444E03EB-B75F-E811-AFBA-F01FAFD8F16A.root
docache /store/mc/RunIIFall17MiniAODv2/WZTo3LNu_3Jets_MLL-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/70000/BC5C9D9E-2965-E811-8F93-E0071B73B6C0.root
