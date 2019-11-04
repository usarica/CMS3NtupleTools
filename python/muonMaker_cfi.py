import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDProducer(
   "MuonMaker",
   aliasprefix      = cms.untracked.string("mus"),

   muonsInputTag    = cms.InputTag("slimmedMuons"),
   vtxInputTag      = cms.InputTag("offlineSlimmedPrimaryVertices"),

   )

