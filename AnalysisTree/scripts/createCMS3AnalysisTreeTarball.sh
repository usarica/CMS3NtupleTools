#!/bin/sh


TARFILE="cms3analysistree.tar"
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
src/MelaAnalytics \
src/ZZMatrixElement/MELA \
src/ZZMatrixElement/setup.sh \
src/CMS3/AnalysisTree \
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

mv $TARFILE $HERE/

popd
