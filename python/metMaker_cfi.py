import FWCore.ParameterSet.Config as cms

metMaker = cms.EDFilter("METMaker",
     metInputTag = cms.InputTag("met"),
     # MC particles
    genParticlesInputTag = cms.InputTag("genParticles")
)


