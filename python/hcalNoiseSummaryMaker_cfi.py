import FWCore.ParameterSet.Config as cms

hcalNoiseSummaryMaker = cms.EDProducer("HcalNoiseSummaryMaker",
	aliasPrefix = cms.untracked.string("hcalnoise"),
                                     hcalNoiseSummaryTag = cms.InputTag("hcalnoise")
)
