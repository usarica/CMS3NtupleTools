import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

muonMaker = cms.EDProducer("MuonMaker",
  aliasPrefix      = cms.untracked.string("mus"),
  muonsInputTag    = cms.InputTag("slimmedMuons"        ),                         
  pfCandsInputTag  = cms.InputTag("packedPFCandidates"),
  vtxInputTag      = cms.InputTag("offlineSlimmedPrimaryVertices"),
  tevMuonsName     = cms.string("tevMuons"),
  cosmicCompat     = cms.InputTag("muons", "cosmicsVeto"),
  pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUp"),
  pfJetsInputTag = cms.InputTag("slimmedJets","",configProcessName.name),
  miniIsoChgValueMap     = cms.InputTag("isoForMu:miniIsoChg"),
  miniIsoAllValueMap     = cms.InputTag("isoForMu:miniIsoAll"),
)


