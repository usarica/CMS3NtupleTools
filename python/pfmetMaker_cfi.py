import FWCore.ParameterSet.Config as cms

pfmetMaker = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETs"),
                            isData              = cms.bool(False)
)


