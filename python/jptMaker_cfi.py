import FWCore.ParameterSet.Config as cms

jptMaker = cms.EDFilter("JPTMaker",
                        jptInputTag       = cms.InputTag("JetPlusTrackZSPCorJetSisCone5"),
                        L2L3jptInputTag   = cms.InputTag("L2L3CorJetSC5JPT"             ),
                        uncorJetsInputTag = cms.InputTag("sisCone5CaloJets"             ),
                        caloJetInputTag   = cms.InputTag("jetMaker", "jetsp4"           )
)
