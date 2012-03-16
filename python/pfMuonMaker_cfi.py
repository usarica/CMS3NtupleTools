import FWCore.ParameterSet.Config as cms

pfMuonMaker = cms.EDProducer(
  "PFMuonMaker",
  pfCandidatesTag = cms.InputTag("particleFlow"),
  isoc_vm_tag     = cms.InputTag("CMS2isoValMuonWithCharged03"),
  ison_vm_tag     = cms.InputTag("CMS2isoValMuonWithNeutral03"),       
  isop_vm_tag     = cms.InputTag("CMS2isoValMuonWithPhotons03"),       
  isoc04_vm_tag   = cms.InputTag("CMS2isoValMuonWithCharged04"),
  ison04_vm_tag   = cms.InputTag("CMS2isoValMuonWithNeutral04"),       
  isop04_vm_tag   = cms.InputTag("CMS2isoValMuonWithPhotons04"),       
  pfAllMuons_tag  = cms.InputTag("pfAllMuons")
)
