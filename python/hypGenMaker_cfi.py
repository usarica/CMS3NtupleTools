import FWCore.ParameterSet.Config as cms

hypGenMaker = cms.EDFilter(
	"HypGenMaker",
	aliasPrefix = cms.untracked.string("hyp"),
	candToGenAssTag = cms.InputTag("candToGenAssMaker"),
	hypInputTag = cms.InputTag("hypDilepMaker")
)

