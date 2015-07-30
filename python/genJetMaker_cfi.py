import FWCore.ParameterSet.Config as cms

genJetMaker = cms.EDProducer("GenJetMaker", 
   genJetsInputTag = cms.InputTag("slimmedGenJets"),
   genJetMinPtCut  = cms.double(10.0),
   aliasPostfix = cms.untracked.string("NoMuNoNu")
)


#genPFJetMaker = cms.EDProducer("GenJetMaker", 
#   genJetsInputTag = cms.InputTag("slimmedGenJets"),
#   genJetMinPtCut  = cms.double(10.0),
#   aliasPostfix = cms.untracked.string("")
#
#)
