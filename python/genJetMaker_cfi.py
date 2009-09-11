import FWCore.ParameterSet.Config as cms

genJetMaker = cms.EDProducer("GenJetMaker", 
   genJetsInputTag = cms.InputTag("sisCone5StGenJets"),
   genJetMinPtCut  = cms.double(10.0)
)
