import FWCore.ParameterSet.Config as cms

sdFilter = cms.EDFilter("SDFilter",
                            elsInputTag_    = cms.InputTag("gsfElectrons"),
                            musInputTag_    = cms.InputTag("muons"),
                            pfjetsInputTag_ = cms.InputTag("ak5PFJets"),
                            photonInputTag_ = cms.InputTag("photons"),
                            metInputTag_    = cms.InputTag("met"),
                            tcmetInputTag_  = cms.InputTag("tcMet"),
                            pfmetInputTag_  = cms.InputTag("pfMet"),
                            elsPt_          = cms.double(10.),   
                            musPt_          = cms.double(10.),   
                            photonPt_       = cms.double(10.),
                            pfjetPt_        = cms.double(40.), 
                            metPt_          = cms.double(60.),   
                            tcmetPt_        = cms.double(60.), 
                            pfmetPt_        = cms.double(60.) 
)
