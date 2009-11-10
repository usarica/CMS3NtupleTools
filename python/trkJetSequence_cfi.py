import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.FastjetParameters_cfi import *
from RecoJets.JetProducers.AntiKtJetParameters_cfi import *

selectTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string('pt > 0.3 & pt<500 & numberOfValidHits > 0')
)
trackCandidates = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("selectTracks"),
    particleType = cms.string('pi+')
)



ak5TrackJets = cms.EDProducer("AntiKtJetProducer",
            FastjetNoPU,
            AntiKtJetParameters,
            CaloJetParameters,
            alias = cms.untracked.string('ANTIKT5TrkJet'),
            FJ_ktRParam = cms.double(0.5)
)
ak5TrackJets.jetPtMin = 1.0
ak5TrackJets.inputEtMin = 0.3
ak5TrackJets.src = cms.InputTag("trackCandidates")
ak5TrackJets.jetType = cms.untracked.string('BasicJet')

cms2TrkJetSequence = cms.Sequence( selectTracks * trackCandidates * ak5TrackJets)
