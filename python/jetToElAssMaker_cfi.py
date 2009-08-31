import FWCore.ParameterSet.Config as cms

jetToElAssMaker = cms.EDFilter("JetToElAssMaker",
                               minDR        = cms.double(0.4)                        ,
                               elInputTag_  = cms.InputTag("electronMaker", "elsp4" ),
                               jetInputTag_ = cms.InputTag("jetMaker"     , "jetsp4")
)


