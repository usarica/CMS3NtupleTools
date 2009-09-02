import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.hltMaker_cfi import *
hltMakerSequence = cms.Sequence(hltMaker8e29*hltMaker1e31)
