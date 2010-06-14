import FWCore.ParameterSet.Config as cms

aSkimFilter = cms.EDProducer("ASkimFilter",
                           filterExpressions = cms.VInputTag(cms.InputTag("electronMaker","elsp4"),
                                                             cms.InputTag("muonMaker","musp4")),
                           filterPtCut  = cms.double(10.0)
)
