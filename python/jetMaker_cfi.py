import FWCore.ParameterSet.Config as cms

jetMaker = cms.EDFilter("JetMaker",
    #EMF corrected Jet Collection
    emfJetCorrectionInputTag = cms.InputTag("L4EMFJetCorJetIcone5"),
    #MC corrected Jet Collection
    mcJetCorrectionInputTag = cms.InputTag("MCJetCorJetIcone5"),
    # jet collection
    jetsInputTag = cms.InputTag("allLayer1Jets"),
    # generator jets collection (for RECO)
    #InputTag genJetsInputTag = midPointCone5GenJets
    # generator jets collection (for AOD)
    genJetsInputTag = cms.InputTag("iterativeCone5GenJets"),
    # MC particles
    genParticlesInputTag = cms.InputTag("genParticles")
)


