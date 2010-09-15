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
                            pfmetPt_        = cms.double(60.), 

                            #L2L3 pfjet correction params
                            doL2L3pfjetCorrection_ = cms.bool(True),
                            PFJetCorrectorL2L3_ = cms.string("ak5PFL2L3"),

                            #thresholds for photon+jet filter
                            photonJet_photonPt_ = cms.double(20.),
                            photonJet_pfjetPt_  = cms.double(30.),
                            photonJet_dr_       = cms.double(0.4)

)
