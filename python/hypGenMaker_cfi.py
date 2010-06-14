import FWCore.ParameterSet.Config as cms

hypGenMaker = cms.EDProducer(
	"HypGenMaker",
	aliasPrefix = cms.untracked.string("hyp"),
	candToGenAssTag = cms.InputTag("candToGenAssMaker"),
	hypInputTag = cms.InputTag("hypDilepMaker")
)

