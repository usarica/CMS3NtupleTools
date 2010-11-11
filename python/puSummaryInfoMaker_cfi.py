import FWCore.ParameterSet.Config as cms

puSummaryInfoMaker = cms.EDProducer("PUSummaryInfoMaker",
                             aliasPrefix = cms.untracked.string("puInfo"),
                             PUInfoInputTag = cms.InputTag("addPileupInfo")
                             )
