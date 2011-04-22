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
                            musPt_          = cms.double(5.),   
                            photonPt_       = cms.double(20.),
                            pfjetPt_        = cms.double(40.), 
                            metPt_          = cms.double(60.),   
                            tcmetPt_        = cms.double(60.), 
                            pfmetPt_        = cms.double(60.),
                            filterName_     = cms.string("nofilter"),
                            tightptcut_        = cms.double(20.),
                            looseptcut_        = cms.double(10.),  

                            #L2L3 pfjet correction params
                            doL2L3pfjetCorrection_ = cms.bool(True),
                            PFJetCorrectorL2L3_ = cms.string("ak5PFL2L3"),

                            #thresholds for photon+jet filter
                            photonJet_photonPt_ = cms.double(20.),
                            photonJet_pfjetPt_  = cms.double(30.),
                            photonJet_dr_       = cms.double(0.4),
                            photonJet_dotrig_   = cms.bool(True),
                        
                            SingleMuTriggerNames_ = cms.untracked.vstring(
                            "HLT_Mu5_v*", 
                            "HLT_Mu8_v*", 
                            "HLT_Mu12_v*",
                            "HLT_Mu8_Jet40_v*",
                            "HLT_Mu30_v*",
                            ),
                            SingleElectronTriggerNames_ = cms.untracked.vstring(
                            "HLT_Ele17_CaloIdL_CaloIsoVL_v*",
                            "HLT_Ele8_v*",
                            "HLT_Ele8_CaloIdL_CaloIsoVL_v*",
                            "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*",
                            "HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*",
                            "HLT_Ele8_CaloIdL_TrkIdVL_v*",
                            ),
                            ElectronHadTriggerNames_ = cms.untracked.vstring(
                            "HLT_DoubleEle8_CaloIdL_TrkIdVL_HT*",
                            "HLT_DoubleEle8_CaloIdT_TrkIdVL_HT*",
                            ),
                            MuHadTriggerNames_ = cms.untracked.vstring(
                            "HLT_DoubleMu3_HT*",
                            "HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT*",
                            "HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT*",
                            ),
                            PhotonTriggerNames_= cms.untracked.vstring(
                             "HLT_Photon75_CaloIdVL_IsoL_v*",
                             "HLT_Photon75_CaloIdVL_v*",
                             "HLT_Photon125_NoSpikeFilter_v*",
                             "HLT_Photon30_CaloIdVL_v*",
                             "HLT_Photon30_CaloIdVL_IsoL_v*",
                             "HLT_Photon20_CaloIdVL_IsoL_v*",
                             "HLT_Photon50_CaloIdVL_IsoL_v*",
                            ),
                            processName = cms.untracked.string("")
)
