import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  pfJetsInputTag                   = cms.InputTag("slimmedJets","","PAT"),
  pfJetPtCut                       = cms.double(5.),
)
ak4JetMaker = cms.EDProducer("AK4JetMaker",
  pfJetsInputTag                   = cms.InputTag("selectedPatJetsAK4PF","","CMS3"),
  pfJetPtCut                       = cms.double(5.),
)
