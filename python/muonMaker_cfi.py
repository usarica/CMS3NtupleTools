import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDProducer("MuonMaker",
  aliasPrefix      = cms.untracked.string("mus"),
  muonsInputTag    = cms.InputTag("slimmedMuons"        ),                         
  beamSpotInputTag = cms.InputTag("beamSpotMaker"),
  pfCandsInputTag  = cms.InputTag("packedPFCandidates"),
  vtxInputTag      = cms.InputTag("offlineSlimmedPrimaryVertices"),
  tevMuonsName     = cms.string("tevMuons"),
  cosmicCompat     = cms.InputTag("muons", "cosmicsVeto"),
  #muonShower       = cms.InputTag("muons", "muonShowerInformation"),
  pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUp")
)


