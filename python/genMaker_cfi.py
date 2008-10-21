import FWCore.ParameterSet.Config as cms

genMaker = cms.EDFilter("GenMaker",
    ntupleOnlyStatus3 = cms.bool(False),
    genParticlesInputTag = cms.InputTag("genParticles")
)


