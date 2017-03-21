import FWCore.ParameterSet.Config as cms

candToGenAssMaker = cms.EDProducer("CandToGenAssExtraMaker",
    #GenJets Input Tag 
    genJetsInputTag = cms.InputTag("slimmedGenJets"),
    #electrons InputTag
    electronsInputTag = cms.InputTag("electronMaker", "elsp4"),
    #phtons Input Tag
    photonsInputTag   = cms.InputTag("photonMaker", "photonsp4"),                                 
    #jets Input Tag
    jetsInputTag = cms.InputTag("jetMaker","jetsp4"),
    #jpts input tag
    pfJetsInputTag = cms.InputTag("pfJetMaker", "pfjetsp4"),                                 
    #ak8jpts input tag
    ak8JetsInputTag = cms.InputTag("subJetMaker", "ak8jetsp4"),                                 
    #muons Input Tag
    muonsInputTag = cms.InputTag("muonMaker","musp4"),
    # MC particles
    genParticlesInputTagPacked = cms.InputTag("packedGenParticles"),
    genParticlesInputTagPruned = cms.InputTag("prunedGenParticles"),
    #jets Input Tag
    tracksInputTag = cms.InputTag("trackMaker", "trkstrkp4"),
    #PIDs of particles to not match to
    #12,14,16->nus, 18->4th gen nu, 1000022->LSP                             
    # must be sorted, since we use binary_search later
    vPIDsToExclude = cms.untracked.vint32(12,14,16,18,1000022)
)


