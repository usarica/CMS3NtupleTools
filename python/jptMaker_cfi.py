import FWCore.ParameterSet.Config as cms

jptMaker = cms.EDFilter("JPTMaker",
                        jptInputTag       = cms.InputTag("JetPlusTrackZSPCorJetIcone5"),
                        L2L3jptInputTag   = cms.InputTag("L2L3CorJetIC5JPT"             ),
                        uncorJetsInputTag = cms.InputTag("prunedUncorrectedCMS2Jets"    )
)
