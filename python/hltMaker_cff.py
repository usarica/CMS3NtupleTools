import FWCore.ParameterSet.Config as cms

from CMS3.NtupleMaker.hltMaker_cfi import *
hltMakerSequence = cms.Sequence(hltMaker)
