import FWCore.ParameterSet.Config as cms

candToGenAssMaker = cms.EDProducer("CandToGenAssMaker",
    #GenJets Input Tag 
    genJetsInputTag = cms.InputTag("cms2antikt5GenJets"),
    #electrons InputTag
    electronsInputTag = cms.InputTag("electronMaker", "elsp4"),
    #phtons Input Tag
    photonsInputTag   = cms.InputTag("photonMaker", "photonsp4"),                                 
    #jets Input Tag
    jetsInputTag = cms.InputTag("jetMaker","jetsp4"),
    #jpts input tag
    pfJetsInputTag = cms.InputTag("pfJetMaker", "pfjetsp4"),                                 
    #muons Input Tag
    muonsInputTag = cms.InputTag("muonMaker","musp4"),
    # MC particles
    genParticlesInputTag = cms.InputTag("genParticles"),
    #jets Input Tag
    tracksInputTag = cms.InputTag("trackMaker", "trkstrkp4"),
    #PIDs of particles to not match to
    #12,14,16->nus, 18->4th gen nu, 1000022->LSP                             
    vPIDsToExclude = cms.untracked.vint32(12,14,16,18,1000022)
)


