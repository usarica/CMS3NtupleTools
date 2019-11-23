#!/bin/sh


TARFILE="cms3ntuplemaker.tar"
echo "SCRAM_ARCH: ${SCRAM_ARCH}"

HERE=$(pwd)

pushd $CMSSW_BASE

tar Jcvf ${TARFILE} \
lib \
biglib \
cfipython \
config \
external \
bin \
src/CMSDataTools \
src/CommonLHETools \
src/ZZMatrixElement \
src/MelaAnalytics \
src/PhysicsTools \
src/RecoEcal \
src/RecoEgamma \
src/EgammaAnalysis \
src/NNKit \
src/CMS3 \
--exclude=lib/${SCRAM_ARCH}/* \
--exclude=src/ZZMatrixElement/MELA/COLLIER/*.so \
--exclude=src/ZZMatrixElement/MELA/data/Pdfdata/NNPDF30_lo_as_0130.LHgrid \
--exclude=src/ZZMatrixElement/MELA/test/reference \
--exclude=src/CMS3/NtupleMaker/*.o \
--exclude=src/CMS3/NtupleMaker/*.so \
--exclude=src/CMS3/NtupleMaker/test/output \
--exclude=src/CMS3/NtupleMaker/test/*.root \
--exclude=src/CMS3/NtupleMaker/test/*.pcm \
--exclude=src/CMS3/NtupleMaker/test/*.so \
--exclude=src/CMS3/NtupleMaker/test/*.a \
--exclude=src/CMS3/NtupleMaker/test/*.o \
--exclude=src/CMS3/NtupleMaker/test/*.d \
--exclude={.git,.gitignore,*.tar,libmcfm*}

mv $TARFILE $HERE/

popd
