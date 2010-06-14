import FWCore.ParameterSet.Config as cms

elToMITConvAssMaker = cms.EDFilter("ElToMITConvAssMaker",
	aliasPrefix = cms.untracked.string("els"),
    elsInputTag = cms.InputTag("electronMaker"),
    mitConvMakerInputTag = cms.InputTag("mitConversionMaker")                               
)
