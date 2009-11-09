import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDFilter("PFJetMaker",
             pfJetsInputTag   = cms.InputTag("sisCone5PFJets"),
             L2L3JetCorrectorName = cms.string("L2L3JetCorrectorSC5PF"),
             pfJetPtCut       = cms.double(5.)  #             
)


