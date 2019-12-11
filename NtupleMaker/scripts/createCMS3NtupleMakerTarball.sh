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
src/ZZMatrixElement/MELA \
src/ZZMatrixElement/setup.sh \
src/MelaAnalytics \
src/PhysicsTools \
src/RecoEcal \
src/RecoEgamma \
src/EgammaAnalysis \
src/NNKit \
src/CMS3 \
--exclude=lib/${SCRAM_ARCH}/* \
--exclude=src/ZZMatrixElement/MELA/data/Pdfdata/NNPDF30_lo_as_0130.LHgrid \
--exclude=src/ZZMatrixElement/MELA/test/reference \
--exclude=src/CommonLHETools/LHEHandler/test/*.root \
--exclude=src/CMS3/AnalysisTree/test/output* \
--exclude=src/CMS3/NtupleMaker/test/output* \
--exclude=src/CMS3/NtupleMaker/test/*.root \
--exclude={.git,.gitignore,__init__.py,*.tar,libmcfm*,*.d,*.a,*.o,*.pcm,*.so,*.pyc,*.mod}

mv $TARFILE $HERE/

popd
