import FWCore.ParameterSet.Config as cms

genJetMaker = cms.EDProducer("GenJetMaker", 
   genJetsInputTag = cms.InputTag("cms2antikt5GenJets"),
   genJetMinPtCut  = cms.double(10.0)
)
