import FWCore.ParameterSet.Config as cms

muonExtraMaker = cms.EDProducer("MuonExtraMaker",
  aliasPrefix      = cms.untracked.string("mus"),
  muonsInputTag    = cms.InputTag("slimmedMuons"        ),                         
  pfCandsInputTag  = cms.InputTag("packedPFCandidates"),
  vtxInputTag      = cms.InputTag("offlineSlimmedPrimaryVertices"),
  tevMuonsName     = cms.string("tevMuons"),
  cosmicCompat     = cms.InputTag("muons", "cosmicsVeto"),
  pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUp")
)


