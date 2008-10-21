import FWCore.ParameterSet.Config as cms

candToGenAssMaker = cms.EDFilter("CandToGenAssMaker",
    #GenJets Input Tag (RECO)
    #InputTag genJetsInputTag =  midPointCone5GenJets
    #GenJets Input Tag (AOD)
    genJetsInputTag = cms.InputTag("iterativeCone5GenJets"),
    #electrons InputTag
    electronsInputTag = cms.InputTag("electronMaker"),
    #jets Input Tag
    jetsInputTag = cms.InputTag("jetMaker"),
    #muons Input Tag
    muonsInputTag = cms.InputTag("muonMaker"),
    # MC particles
    genParticlesInputTag = cms.InputTag("genParticles"),
    #jets Input Tag
    tracksInputTag = cms.InputTag("trackMaker")
)


