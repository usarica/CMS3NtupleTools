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
cfipython \
config \
external \
bin \
src/CMSDataTools \
src/MelaAnalytics \
src/ZZMatrixElement/MELA \
src/ZZMatrixElement/setup.sh \
src/CMS3/AnalysisTree \
src/CMS3/Dictionaries \
src/CMS3/MELAHelpers \
--exclude=lib/${SCRAM_ARCH}/* \
--exclude=src/ZZMatrixElement/MELA/test/reference \
--exclude=src/CMS3/AnalysisTree/test/output* \
--exclude=src/CMS3/AnalysisTree/test/manual* \
--exclude=src/CMS3/AnalysisTree/test/Pdfdata \
--exclude=src/CMS3/AnalysisTree/test/*.root \
--exclude=src/CMS3/AnalysisTree/test/br.sm* \
--exclude=src/CMS3/AnalysisTree/test/*.dat \
--exclude=src/CMS3/AnalysisTree/test/*.DAT \
--exclude=src/CMS3/AnalysisTree/test/*.txt \
--exclude=src/CMS3/AnalysisTree/test/*.csv \
--exclude=src/CMS3/AnalysisTree/test/*.md \
--exclude=src/CMS3/AnalysisTree/test/*.sh \
--exclude={.git,.gitignore,__init__.py,*.tar,libmcfm*,*.d,*.a,*.o,*.pcm,*.so,*.pyc,*.mod}


if [[ ! -z "$extraTarFile" ]];then
  rm -f $extraTarFile
fi

mv $TARFILE $HERE/

popd
