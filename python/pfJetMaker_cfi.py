import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  aliasPrefix = cms.untracked.string("pfjets"),
  pfJetsInputTag                   = cms.InputTag("slimmedJets","","PAT"),
  pfJetPtCut                       = cms.double(5.),
)
pfJetPUPPIMaker = cms.EDProducer("PFJetMaker",
  aliasPrefix = cms.untracked.string("pfjets_puppi"),
  pfJetsInputTag                   = cms.InputTag("slimmedJetsPuppi","","PAT"),
  pfJetPtCut                       = cms.double(20.), #miniAOD doesn't go lower than 20
)
ak4JetMaker = cms.EDProducer("AK4JetMaker",
  pfJetsInputTag                   = cms.InputTag("selectedPatJetsAK4PF","","CMS3"),
  pfJetPtCut                       = cms.double(5.),
)
