import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.SISConeJetParameters_cfi import *
from RecoJets.JetProducers.FastjetParameters_cfi import *

selectTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string('pt > 0.3 & pt<500 & numberOfValidHits > 0')
)
trackCandidates = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("selectTracks"),
    particleType = cms.string('pi+')
)
scTrackJetParameters = cms.PSet(
    SISConeJetParameters,FastjetNoPU, 
    jetType = cms.untracked.string('BasicJet'),
    src = cms.InputTag("trackCandidates"),
    jetPtMin = cms.double( 0.0 ),
    inputEtMin = cms.double(0.3),
    inputEMin = cms.double(0.0)
)
SISCone5TrkJets = cms.EDProducer("SISConeJetProducer",
    scTrackJetParameters,
    alias = cms.untracked.string('SISCone5TrkJets'),
    coneRadius = cms.double(0.5)
)


cms2TrkJetSequence = cms.Sequence( selectTracks * trackCandidates * SISCone5TrkJets)
