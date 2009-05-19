import FWCore.ParameterSet.Config as cms

pfmetMaker = cms.EDFilter("PFMETMaker",
                          pfMetInputTag_ = cms.InputTag("pfMet")
)


