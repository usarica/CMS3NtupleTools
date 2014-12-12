import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.electronMaker_cfi import *
from CMS2.NtupleMaker.muonMaker_cfi import *

aSkimFilter = cms.EDFilter("ASkimFilter",
                           electronsInputTag = electronMaker.electronsInputTag,
                           muonsInputTag     = muonMaker.muonsInputTag,
                           useSTAMuons       = cms.bool(False),
                           eleFilterPtCut    = cms.double(10.0),
                           muFilterPtCut     = cms.double(5.0)
)
