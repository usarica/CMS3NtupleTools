import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

hcalNoiseSummaryMaker = cms.EDProducer("HcalNoiseSummaryMaker",
	aliasPrefix = cms.untracked.string("hcalnoise"),
                                     hcalNoiseSummaryTag = cms.InputTag("hcalnoise","",configProcessName.name)
)
