import FWCore.ParameterSet.Config as cms

hcalNoiseSummaryMaker = cms.EDFilter("HcalNoiseSummaryMaker",
                                     hcalNoiseSummaryTag = cms.InputTag("hcalnoise")
)
