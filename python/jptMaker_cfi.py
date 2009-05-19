import FWCore.ParameterSet.Config as cms

jptMaker = cms.EDFilter("JPTMaker",
                        jptInputTag      = cms.InputTag("JetPlusTrackZSPCorJetSisCone5"),
                        L2L3jptInputTag  = cms.InputTag("L2L3CorJetIC5JPT"             ),
                        caloJetsInputTag = cms.InputTag("L2L3CorJet"                   )
)
