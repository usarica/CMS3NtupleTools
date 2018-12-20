
import FWCore.ParameterSet.Config as cms

pfmetMaker = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETs"),
                            onlySaveTwoVector   = cms.bool(False),
                            doUncertainties   = cms.bool(True)
)

pfmetpuppiMaker = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt_puppi"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETsPuppi"),
                            onlySaveTwoVector   = cms.bool(False),
                            doUncertainties   = cms.bool(True)
)


