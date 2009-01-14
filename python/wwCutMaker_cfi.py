import FWCore.ParameterSet.Config as cms

wwCutMaker = cms.EDFilter("WWCutMaker",
       dileptonInputTag = cms.InputTag("hypDilepMaker"),
       muonsInputTag = cms.InputTag("muonMaker"),
       electronsInputTag = cms.InputTag("electronMaker"),
       eltomuassInputTag = cms.InputTag("elToMuAssMaker"),
       genpsInputTag = cms.InputTag("genMaker")
)
