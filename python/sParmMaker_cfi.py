import FWCore.ParameterSet.Config as cms

sParmMaker = cms.EDProducer("SParmMaker",
                            sparm_inputTag       = cms.InputTag("source"),
                            aliasPrefix          = cms.untracked.string("sparm"),
                            vsparms               = cms.untracked.vstring()
)
