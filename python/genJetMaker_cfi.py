import FWCore.ParameterSet.Config as cms

genJetMaker = cms.EDProducer("GenJetMaker", 
   genJetsInputTag = cms.InputTag("cms2antikt5GenJets"),
   genJetMinPtCut  = cms.double(10.0),
   aliasPostfix = cms.untracked.string("")
)


genPFJetMaker = cms.EDProducer("GenJetMaker", 
   genJetsInputTag = cms.InputTag("cms2antikt5PFGenJets"),
   genJetMinPtCut  = cms.double(10.0),
   aliasPostfix = cms.untracked.string("NoNu")

)
