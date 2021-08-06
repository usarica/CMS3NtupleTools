#!/bin/bash

(
set -euo pipefail

cd $(dirname ${BASH_SOURCE[0]})

strCMSSW_VERSION="$(echo ${CMSSW_VERSION} | sed -e s/CMSSW_// -e s/_//g -e s/patch\.//)"
strCMSSW_VERSION_MAJOR="$(echo ${CMSSW_VERSION} | sed -e s/CMSSW_// -e s/_.*//g)"
strCMSSW_VERSION_MINOR="$(echo ${CMSSW_VERSION} | sed -e s/CMSSW_// -e 's/_/ /g' | awk '{print $$2}')"
echo "CMSSW version: ${strCMSSW_VERSION}"
echo "CMSSW major version: ${strCMSSW_VERSION_MAJOR}"
echo "CMSSW minor version: ${strCMSSW_VERSION_MINOR}"

declare -a setupArgs=()

for farg in "$@"; do
  setupArgs+=( "$farg" ) 
done
declare -i nSetupArgs
nSetupArgs=${#setupArgs[@]}
if [[ $nSetupArgs -eq 0 ]]; then
  setupArgs+=( -j 1 )
  nSetupArgs=2
fi

if [[ "$nSetupArgs" -eq 1 ]] && [[ "${setupArgs[0]}" == *"clean"* ]]; then
  rm -f BuildFile.xml
  rm -f plugins/BuildFile.xml
elif [[ "$nSetupArgs" -ge 1 ]] && [[ "$nSetupArgs" -le 2 ]] && [[ "${setupArgs[0]}" == *"-j"* ]]; then
  rm -f BuildFile.xml
  rm -f plugins/BuildFile.xml
  cp BuildFile.tpl BuildFile.xml
  cp plugins/BuildFile.tpl plugins/BuildFile.xml
  sed -i "s|.oOCMSSW_VERSIONOo.|${strCMSSW_VERSION}|g" BuildFile.xml
  sed -i "s|.oOCMSSW_VERSION_MAJOROo.|${strCMSSW_VERSION_MAJOR}|g" BuildFile.xml
  sed -i "s|.oOCMSSW_VERSION_MINOROo.|${strCMSSW_VERSION_MINOR}|g" BuildFile.xml
  sed -i "s|.oOCMSSW_VERSIONOo.|${strCMSSW_VERSION}|g" plugins/BuildFile.xml
  sed -i "s|.oOCMSSW_VERSION_MAJOROo.|${strCMSSW_VERSION_MAJOR}|g" plugins/BuildFile.xml
  sed -i "s|.oOCMSSW_VERSION_MINOROo.|${strCMSSW_VERSION_MINOR}|g" plugins/BuildFile.xml
else
  echo "Unknown arguments:"
  echo "  ${setupArgs[@]}"
  exit 1
fi

scramv1 b "${setupArgs[@]}"

)
