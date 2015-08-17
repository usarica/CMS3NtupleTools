import FWCore.ParameterSet.Config as cms

METToolboxJetMaker = cms.EDProducer("T1CHSMETJets",
  aliasPrefix     = cms.untracked.string("pfjets_METToolbox"),
  pfJetsInputTag  = cms.InputTag("ak4PFJetsCHS","","CMS3"),
  pfCandidatesTag = cms.InputTag("packedPFCandidates"),
  pfJetPtCut      = cms.double(0.),
)
