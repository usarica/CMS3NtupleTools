import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5TrackJets_cfi import ak5TrackJets
from RecoJets.JetProducers.TracksForJets_cff import *

cms2TrkJetSequence = cms.Sequence(trackWithVertexRefSelector*trackRefsForJets*ak5TrackJets)
