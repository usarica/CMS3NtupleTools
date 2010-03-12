import FWCore.ParameterSet.Config as cms

hcalNoiseSummaryMaker = cms.EDFilter("HcalNoiseSummaryMaker",
	aliasPrefix = cms.untracked.string("hcalnoise"),
                                     hcalNoiseSummaryTag = cms.InputTag("hcalnoise")
)
