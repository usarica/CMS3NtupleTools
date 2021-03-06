import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDProducer(
   "MuonMaker",
   aliasprefix = cms.untracked.string("muons"),

   year = cms.int32(-1), # Must be overriden by main_pset

   refurbishSelections = cms.bool(False),

   muonsInputTag = cms.InputTag("slimmedMuons"),

   pfCandidatesInputTag = cms.InputTag("packedPFCandidates"),

   vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

   rhoInputTag = cms.InputTag("fixedGridRhoFastjetAll"),

   )

