import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  pfJetsInputTag                   = cms.InputTag("slimmedJets","","PAT"),
  #pfCandidatesTag                  = cms.InputTag("packedPFCandidates","","PAT"),
  pfJetPtCut                       = cms.double(5.),
  #PFJetCorrectorL2L3               = cms.InputTag("ak4PFCHSL2L3"),
  #PFJetCorrectorL1Fast             = cms.InputTag("ak4PFCHSL1Fastjet"),
  #PFJetCorrectorL1FastL2L3         = cms.InputTag("ak4PFCHSL1FastL2L3")
  PFJetCorrectorL2L3               = cms.string("ak4PFCHSL2L3"),
  PFJetCorrectorL1Fast             = cms.string("ak4PFCHSL1Fastjet"),
  PFJetCorrectorL1FastL2L3         = cms.string("ak4PFCHSL1FastL2L3")
)
