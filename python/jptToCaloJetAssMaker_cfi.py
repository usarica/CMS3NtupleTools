import FWCore.ParameterSet.Config as cms


jptToCaloJetAssMaker = cms.EDFilter("JPTtoCaloJetAssMaker",
                               cms2CaloJetInputTag  = cms.InputTag("jetMaker", "jetsp4" ),
                               jptInputTag = cms.InputTag("jptMaker", "jptsp4")
)
