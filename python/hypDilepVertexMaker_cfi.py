import FWCore.ParameterSet.Config as cms

hypDilepVertexMaker = cms.EDFilter("HypDilepVertexMaker",
	aliasPrefix = cms.untracked.string("hyp"),
                                   recomuonsInputTag = cms.InputTag("muons"),
                                   cms2muonsInputTag = cms.InputTag("muonMaker"),
                                   recoelectronsInputTag = cms.InputTag("gsfElectrons"),
                                   cms2electronsInputTag = cms.InputTag("electronMaker"),
                                   hypInputTag = cms.InputTag("hypDilepMaker")
                                   )
