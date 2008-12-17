import FWCore.ParameterSet.Config as cms

jetMaker = cms.EDFilter("JetMaker",
    #EMF corrected Jet Collection
    emfJetCorrectionInputTag = cms.InputTag("L2L3L4CorJet"),
    #MC corrected Jet Collection
    mcJetCorrectionInputTag = cms.InputTag("L2L3CorJet"),
    # jet collection
    jetsInputTag = cms.InputTag("selectedLayer1Jets"),
    # generator jets collection (for RECO)
    #genJetsInputTag = cms.InputTag("sisCone5GenJets"), 
    # generator jets collection (for AOD)
    #genJetsInputTag = cms.InputTag("sisCone5GenJets"),
    # MC particles
    #genParticlesInputTag = cms.InputTag("genParticles")
)


