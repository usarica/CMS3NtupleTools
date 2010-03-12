import FWCore.ParameterSet.Config as cms

pfmetMaker = cms.EDFilter("PFMETMaker",
	aliasPrefix = cms.untracked.string("evt"),
                          pfMetInputTag_ = cms.InputTag("pfMet")
)


