import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  pfJetsInputTag                   = cms.InputTag("slimmedJets"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(5.),
  PFJetCorrectorL2L3               = cms.string("ak4PFCHSL2L3"),
  PFJetCorrectorL1Fast             = cms.string("ak4PFCHSL1Fastjet"),
  PFJetCorrectorL1FastL2L3         = cms.string("ak4PFCHSL1FastL2L3")
)
