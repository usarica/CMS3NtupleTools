import FWCore.ParameterSet.Config as cms

luminosityMaker = cms.EDProducer("LuminosityMaker",
                                 aliasPrefix = cms.untracked.string("ls"),
                                 lumiSummaryInputTag = cms.InputTag("lumiProducer")
                                 )


