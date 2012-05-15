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
                            photonJet_pfjetPt_hzz_ = cms.double(50.),
                            photonJet_doJet_hzz_ = cms.bool(True),
                            photonJet_dr_       = cms.double(0.4),
                            photonJet_dotrig_   = cms.bool(True),
                        
                            SingleMuTriggerNames_ = cms.untracked.vstring(
                            "HLT_Mu*", 
                            "HLT_IsoMu*", 
                            "HLT_RelIso1p0*", 
                            ),
                            SingleElectronTriggerNames_ = cms.untracked.vstring(
                            "HLT_Ele25*", #temporary --> may remove/restrict later
                            ),
                            DoubleElectronTriggerNames_ = cms.untracked.vstring(
                            "HLT_Ele17*",
                            "HLT_Photon*",
                            "HLT_Ele8*",
                            "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*",
                            "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*",
                            ),
                            ElectronHadTriggerNames_ = cms.untracked.vstring(
                            "HLT_DoubleEle*",
                            "HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*",
                            "HLT_Ele25*", #temporary --> may remove/restrict later
                            "HLT_Ele27*", #temporary --> may remove/restrict later
                            ),
                            MuHadTriggerNames_ = cms.untracked.vstring(
                            "HLT*Mu*Ele*HT*",
                            "HLT*Doubl*Mu*HT*",
                            "HLT_IsoMu20*", #temporary --> may remove/restrict later
                            "HLT_Mu20*",    #temporary --> may remove/restrict later
                            ),
                            PhotonTriggerNames_= cms.untracked.vstring(
                             "HLT_Photon20_CaloIdVL_IsoL_v*",
                             #"HLT_Photon30_CaloIdVL_v*",
                             "HLT_Photon30_CaloIdVL_IsoL_v*",
                             "HLT_Photon50_CaloIdVL_IsoL_v*",
                             "HLT_Photon75_CaloIdVL_IsoL_v*",
                             #"HLT_Photon75_CaloIdVL_v*",
                             "HLT_Photon90_CaloIdVL_IsoL_v*",
                             #"HLT_Photon125_v*",
                             #"HLT_Photon125_NoSpikeFilter_v*",
                             #"HLT_Photon135_v*",
                             #"HLT_Photon400_v*",
                             "HLT_Photon135_v*",
                             "HLT_Photon150_v*",
                             "HLT_Photon160_v*",
                            ),
                            processName = cms.untracked.string("")
)
