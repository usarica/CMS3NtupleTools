import FWCore.ParameterSet.Config as cms

jptMaker = cms.EDFilter("JPTMaker",
                        jptInputTag       = cms.InputTag("JetPlusTrackZSPCorJetIcone5")
)
