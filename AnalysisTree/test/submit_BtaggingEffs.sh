#!/bin/bash

date=$1
period=$2
prodVersion=$3
function=produceBtaggingEfficiencies
script=${function}.cc
jobdate="BtaggingEfficiencies_${date}"
arguments='"<strSampleSet>","<period>","<prodVersion>","<strdate>",<ichunk>,<nchunks>,<theGlobalSyst>,<applyTightLeptonVetoIdToAK4Jets>'
arguments="${arguments/<period>/$period}"
arguments="${arguments/<prodVersion>/$prodVersion}"
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<applyTightLeptonVetoIdToAK4Jets>/false}"

skimdir="/ceph/cms/store/user/usarica/Offshell_2L2Nu/Skims/${prodVersion}"
skimdir=$( ExecuteCompiledCommand GetStandardHostPathToStore ${skimdir} t2.ucsd.edu )
if [[ "$( ExecuteCompiledCommand DirectoryExists ${skimdir} )" == "false" ]]; then
  exit 1
fi

MCSysts=( sNominal eJECDn eJECUp eJERDn eJERUp ePUDn ePUUp ePUJetIdEffDn ePUJetIdEffUp )

for d in $(ls $skimdir); do
  for dd in $(ls $skimdir/$d); do
    if [[ $dd == "Run20"* ]];then
      continue
    fi
    sample="/$d/$dd"
    sampledir=${skimdir}${sample}
    sample="${sample}/MINIAODSIM"

    let nfiles=$(ExecuteCompiledCommand lsdir $sampledir | grep -e ".root" | wc -l)
    echo "$sampledir has $nfiles files"

    nchunks=$((nfiles / 4))
    nfeff=$((nchunks * 4))
    if [[ $nfeff -lt $nfiles ]];then
      nchunks=$((nchunks + 1))
    fi
    echo " - There will be $nchunks chunks"
    if [[ $nchunks -eq 0 ]]; then
      echo "Oh-uh..."
    fi

    for syst in "${MCSysts[@]}"; do

      let ichunk=0
      while [[ $ichunk -lt $nchunks ]]; do
        strargs="${arguments}"
        strargs="${strargs/<strSampleSet>/${sample}}"
        strargs="${strargs/<ichunk>/${ichunk}}"
        strargs="${strargs/<nchunks>/${nchunks}}"
        strargs="${strargs/<theGlobalSyst>/${syst}}"


        submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"

        let ichunk=$ichunk+1
      done

    done
  done
  
  echo "====="
done
