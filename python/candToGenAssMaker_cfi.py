import FWCore.ParameterSet.Config as cms

candToGenAssMaker = cms.EDFilter("CandToGenAssMaker",
    #GenJets Input Tag 
    genJetsInputTag = cms.InputTag("antikt5GenJets"),
    #electrons InputTag
    electronsInputTag = cms.InputTag("electronMaker"),
    #jets Input Tag
    jetsInputTag = cms.InputTag("jetMaker"),
    #muons Input Tag
    muonsInputTag = cms.InputTag("muonMaker"),
    # MC particles
    genParticlesInputTag = cms.InputTag("genParticles"),
    #jets Input Tag
    tracksInputTag = cms.InputTag("trackMaker"),
    #PIDs of particles to not match to
    #12,14,16->nus, 18->4th gen nu, 1000022->LSP                             
    vPIDsToExclude = cms.untracked.vint32(12,14,16,18,1000022)
)


