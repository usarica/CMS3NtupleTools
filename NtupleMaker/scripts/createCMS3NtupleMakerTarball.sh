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
src/JHUGenMELA/MELA \
src/JHUGenMELA/setup.sh \
src/MelaAnalytics \
src/PhysicsTools \
src/RecoEcal \
src/RecoEgamma \
src/EgammaAnalysis \
src/NNKit \
src/CMS3 \
--exclude=lib/${SCRAM_ARCH}/* \
--exclude=src/JHUGenMELA/MELA/test/reference \
--exclude=src/CommonLHETools/LHEHandler/test/*.root \
--exclude=src/CMS3/AnalysisTree \
--exclude=src/CMS3/NtupleMaker/test/output* \
--exclude=src/CMS3/NtupleMaker/test/manual* \
--exclude=src/CMS3/NtupleMaker/test/Pdfdata \
--exclude=src/CMS3/NtupleMaker/test/*.root \
--exclude=src/CMS3/NtupleMaker/test/br.sm* \
--exclude=src/CMS3/NtupleMaker/test/*.dat \
--exclude=src/CMS3/NtupleMaker/test/*.DAT \
--exclude=src/CMS3/NtupleMaker/test/*.txt \
--exclude=src/CMS3/NtupleMaker/test/*.csv \
--exclude=src/CMS3/NtupleMaker/test/*.md \
--exclude=src/CMS3/NtupleMaker/test/*.sh \
--exclude={.git,.gitignore,__init__.py,*.tar,libmcfm*,*.d,*.a,*.o,*.pcm,*.so,*.pyc,*.mod,*.bkp}

mv $TARFILE $HERE/

popd
