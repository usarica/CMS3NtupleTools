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

tar Jcvf ${TARFILE} ${extraTarFile} \
lib \
biglib \
bin \
src/CMSDataTools \
src/MelaAnalytics \
src/ZZMatrixElement/MELA \
src/ZZMatrixElement/setup.sh \
src/CMS3/AnalysisTree \
src/CMS3/Dictionaries \
src/CMS3/MELAHelpers \
--exclude=src/*/*/src \
--exclude=src/*/*/bin \
--exclude=src/*/*/scripts \
--exclude=src/ZZMatrixElement/MELA/test/reference \
--exclude=src/ZZMatrixElement/MELA/test/Pdfdata \
--exclude=src/ZZMatrixElement/MELA/test/*.d \
--exclude=src/ZZMatrixElement/MELA/test/*.a \
--exclude=src/ZZMatrixElement/MELA/test/*.o \
--exclude=src/ZZMatrixElement/MELA/test/*.so \
--exclude=src/ZZMatrixElement/MELA/test/*.pcm \
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
--exclude={.git,.gitignore,__init__.py,*.tar,*.pyc,*.mod,*.out}


if [[ ! -z "$extraTarFile" ]];then
  rm -f $extraTarFile
fi

mv $TARFILE $HERE/

popd
