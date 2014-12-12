import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.trkMetMaker_cfi import *

myTrkMetMaker = trkMetMaker.clone()

myTrkJetMetMaker = trkMetMaker.clone()
myTrkJetMetMaker.correctJet = cms.bool(True)
myTrkJetMetMaker.aliasPrefix = cms.untracked.string("trkjet")

trkMetSequence = cms.Sequence( myTrkMetMaker + myTrkJetMetMaker )

