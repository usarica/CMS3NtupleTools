import FWCore.ParameterSet.Config as cms

pfmetMaker = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt"),
                            pfMetInputTag_ = cms.InputTag("pfMet"),
                            pfMetCorInputTag_ = cms.InputTag("pfType1CorrectedMet")
)


