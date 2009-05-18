import FWCore.ParameterSet.Config as cms

jptMaker = cms.EDFilter("JPTMaker",
                        jptInputTag      = cms.InputTag("JetPlusTrackZSPCorJetSisCone"),
                        caloJetsInputTag = cms.InputTag("sisCone5CaloJets"            )
)


