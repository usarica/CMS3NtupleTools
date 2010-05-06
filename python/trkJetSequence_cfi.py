import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.trackJetCollectionPruner_cfi import *
from CMS2.NtupleMaker.trkJetMaker_cfi import *

cms2TrkJetSequence = cms.Sequence(prunedUncorrectedCMS2TrackJets * trkJetMaker)
