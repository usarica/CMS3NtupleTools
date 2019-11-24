import FWCore.ParameterSet.Config as cms

sParmMaker = cms.EDProducer("SParmMaker",
                            sparm_inputTag       = cms.InputTag("externalLHEProducer"),
                            config_inputTag       = cms.InputTag("generator"),
                            aliasPrefix          = cms.untracked.string("sparm"),
                            vsparms               = cms.untracked.vstring()
)
