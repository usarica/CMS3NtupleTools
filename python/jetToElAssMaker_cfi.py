import FWCore.ParameterSet.Config as cms

jetToElAssMaker = cms.EDFilter("JetToElAssMaker",
                               minDR       = cms.double(0.1)                        ,
                               elInputTag  = cms.InputTag("electronMaker", "elsp4" ),
                               jetInputTag = cms.InputTag("jetMaker"     , "jetsp4")
)


