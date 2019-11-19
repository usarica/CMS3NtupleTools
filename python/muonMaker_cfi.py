import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDProducer(
   "MuonMaker",
   aliasprefix      = cms.untracked.string("mus"),
   year = cms.int32(-1), # Must be overriden by main_pset

   muonsInputTag    = cms.InputTag("slimmedMuons"),
   vtxInputTag      = cms.InputTag("offlineSlimmedPrimaryVertices"),

   rhoInputTag = cms.InputTag("fixedGridRhoFastjetAll"),

   )

