#!/bin/sh


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


TARFILE="cms3analysistree.tar"
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

extraIncludes=""
extraExcludes=""
if [[ -d ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit ]]; then
  extraIncludes="${extraIncludes} src/HiggsAnalysis/CombinedLimit"
  extraExcludes="${extraExcludes} --exclude=src/HiggsAnalysis/CombinedLimit/data --exclude=src/HiggsAnalysis/CombinedLimit/doc*"

  echo "extraIncludes is now ${extraIncludes}"
  echo "extraExcludes is now ${extraExcludes}"
fi

tar Jcvf ${TARFILE} ${extraTarFile} ${extraIncludes} \
lib \
biglib \
bin \
src/IvyFramework \
src/MelaAnalytics \
src/JHUGenMELA/MELA \
src/JHUGenMELA/setup.sh \
src/CMS3/AnalysisTree \
src/CMS3/Dictionaries \
--exclude=src/*/*/src \
--exclude=src/*/*/bin \
--exclude=src/*/*/scripts ${extraExcludes} \
--exclude=src/JHUGenMELA/MELA/test/reference \
--exclude=src/JHUGenMELA/MELA/test/Pdfdata \
--exclude=src/JHUGenMELA/MELA/test/*.d \
--exclude=src/JHUGenMELA/MELA/test/*.a \
--exclude=src/JHUGenMELA/MELA/test/*.o \
--exclude=src/JHUGenMELA/MELA/test/*.so \
--exclude=src/JHUGenMELA/MELA/test/*.pcm \
--exclude=src/*/*/test/Pdfdata \
--exclude=src/*/*/test/br.sm* \
--exclude=src/*/*/test/*.dat \
--exclude=src/*/*/test/*.DAT \
--exclude=src/*/*/test/tmp* \
--exclude=src/*/*/test/temp* \
--exclude=src/CMS3/AnalysisTree/test/output* \
--exclude=src/CMS3/AnalysisTree/test/manual* \
--exclude=src/CMS3/AnalysisTree/test/*.root \
--exclude=src/CMS3/AnalysisTree/test/*.txt \
--exclude=src/CMS3/AnalysisTree/test/*.csv \
--exclude=src/CMS3/AnalysisTree/test/*.md \
--exclude=src/CMS3/AnalysisTree/test/*.sh \
--exclude=src/CMS3/AnalysisTree/test/*.d \
--exclude=src/CMS3/AnalysisTree/test/*.a \
--exclude=src/CMS3/AnalysisTree/test/*.o \
--exclude=src/CMS3/AnalysisTree/test/*.so \
--exclude=src/CMS3/AnalysisTree/test/*.pcm \
--exclude={.git,.gitignore,__init__.py,*.tar,*.pyc,*.mod,*.out,*.bkp}


if [[ ! -z "$extraTarFile" ]];then
  rm -f $extraTarFile
fi

mv $TARFILE $HERE/

popd
