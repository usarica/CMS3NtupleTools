#!/bin/sh


if [[ ! -d ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit ]]; then
  echo "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit must be present for this tarball."
  exit 1
fi

EXTRAFILES=""
for fargo in "$@";do
  fcnargname=""
  farg="${fargo//\"}"
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "addfile="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    fcnargname="$(readlink -e $fcnargname)"
    removename="${CMSSW_BASE}/"
    fcnargname="${fcnargname//$removename}"
    if [[ -z "$EXTRAFILES" ]];then
      EXTRAFILES="$fcnargname"
    else
      EXTRAFILES="$EXTRAFILES $fcnargname"
    fi
  fi
done


TARFILE="cms3tnpfit.tar"
echo "SCRAM_ARCH: ${SCRAM_ARCH}"

HERE=$(pwd)

pushd $CMSSW_BASE

extraTarFile=""
if [[ ! -z "$EXTRAFILES" ]];then
  extraTarFile="extras.tar"
  echo "Will include these files in ${extraTarFile}: ${EXTRAFILES}"
  tar Jcvf $extraTarFile $EXTRAFILES
  mv $extraTarFile extras.more # HACK: Hide file extension
  extraTarFile="extras.more"
fi

tar Jcvf ${TARFILE} ${extraTarFile} \
lib/*/*CombinedLimit* \
lib/*/libCMSDataToolsAnalysisTree.so \
lib/*/libJHUGenMELAMELA.so \
biglib \
bin \
src/HiggsAnalysis/CombinedLimit \
src/JHUGenMELA/MELA/data/${SCRAM_ARCH}/lib*.so \
--exclude=src/*/*/src \
--exclude=src/*/*/bin \
--exclude=src/*/*/scripts \
--exclude=src/HiggsAnalysis/CombinedLimit/data \
--exclude=src/HiggsAnalysis/CombinedLimit/macros \
--exclude=src/HiggsAnalysis/CombinedLimit/doc* \
--exclude=src/*/*/test/Pdfdata \
--exclude=src/*/*/test/br.sm* \
--exclude=src/*/*/test/*.dat \
--exclude=src/*/*/test/*.DAT \
--exclude=src/*/*/test/tmp* \
--exclude=src/*/*/test/temp* \
--exclude={.git,.gitignore,__init__.py,*.tar,*.pyc,*.mod,*.out}


if [[ ! -z "$extraTarFile" ]];then
  rm -f $extraTarFile
fi

mv $TARFILE $HERE/

popd
