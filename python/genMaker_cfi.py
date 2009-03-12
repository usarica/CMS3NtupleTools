import FWCore.ParameterSet.Config as cms

genMaker = cms.EDFilter(
	"GenMaker",
    ntupleOnlyStatus3 = cms.bool(True),
	ntupleDaughters = cms.bool(True),
    genParticlesInputTag = cms.InputTag("genParticles")
)


