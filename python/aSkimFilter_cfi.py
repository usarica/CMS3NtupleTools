import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.electronMaker_cfi import *
from CMS2.NtupleMaker.muonMaker_cfi import *

aSkimFilter = cms.EDFilter("ASkimFilter",
                             #electronsInputTag = cms.InputTag("gsfElectrons"),
                             electronsInputTag = electronMaker.electronsInputTag,
                             #muonsInputTag     = cms.InputTag("muons"        ),
                             muonsInputTag     = muonMaker.muonsInputTag,
                             useSTAMuons       = cms.bool(False),
                             filterPtCut       = cms.double(10.0)
)
