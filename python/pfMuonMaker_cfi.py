import FWCore.ParameterSet.Config as cms

pfMuonMaker = cms.EDProducer("PFMuonMaker",
                             pfCandidatesTag = cms.InputTag("particleFlow"),
                             isoc_vm_tag     = cms.InputTag("CMS2isoValMuonWithCharged"),
                             ison_vm_tag     = cms.InputTag("CMS2isoValMuonWithNeutral"),       
                             isop_vm_tag     = cms.InputTag("CMS2isoValMuonWithPhotons"),       
                             isoc04_vm_tag   = cms.InputTag("CMS2isoValMuonWithCharged"),
                             ison04_vm_tag   = cms.InputTag("CMS2isoValMuonWithNeutral"),       
                             isop04_vm_tag   = cms.InputTag("CMS2isoValMuonWithPhotons"),       
                             pfAllMuons_tag  = cms.InputTag("pfAllMuons")
)
