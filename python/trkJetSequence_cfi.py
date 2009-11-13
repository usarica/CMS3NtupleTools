import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.TracksForJets_cff import *
from RecoJets.JetProducers.ak5TrackJets_cfi import *

cms2TrkJetSequence = cms.Sequence( tracksForJets * ak5TrackJets)
