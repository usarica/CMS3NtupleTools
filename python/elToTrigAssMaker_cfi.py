import FWCore.ParameterSet.Config as cms

elToTrigAssMaker = cms.EDProducer("ObjectToTriggerLegAssMaker",

    aliasPrefix         = cms.untracked.string("els"),
    objectInputTag      = cms.untracked.InputTag("electronMaker:elsp4"),
    cone                = cms.untracked.double(0.2),

    triggers    = cms.untracked.VInputTag(
        ### RUN II ###
        # Double-lepton matching and single-lepton matching can be done using the TriggerObjects in CMS3, for example, 
        # CORE::TriggerSelections::passUnprescaledHLTTrigger(const char* arg, const LorentzVector &obj)
        # But for LeadingLeg matching we need to know the specific filter passed by the leading leg
        # Use:    'TriggerName:FilterName:CMS3name'   --> this requires that electron matches to trigger object passing FilterName of TriggerName
        # Use:    'TriggerName::CMS3name'             --> this requires that electron matches to trigger object passing the LAST EDFILTER of TriggerName
        # cms.InputTag('HLT_Ele25WP60_Ele8_Mass55_v*:hltEle25WP60Ele8TrackIsoFilter:HLT_Ele25WP60_Ele8_Mass55_LeadingLeg'),
        # cms.InputTag('HLT_Ele25WP60_Ele8_Mass55_v*::HLT_Ele25WP60_Ele8_Mass55'),

        # cms.InputTag('HLT_Ele25WP60_SC4_Mass55_v*:hltEle25WP60SC4TrackIsoFilter:HLT_Ele25WP60_SC4_Mass55_LeadingLeg'),
        # cms.InputTag('HLT_Ele25WP60_SC4_Mass55_v*::HLT_Ele25WP60_SC4_Mass55'),

        # cms.InputTag('HLT_Ele5_SC5_JPsi_Mass2to4p5_v*:hltEle5SC5JPsiTrackIsoFilter:HLT_Ele5_SC5_JPsi_Mass2to4p5_LeadingLeg'),
        # cms.InputTag('HLT_Ele5_SC5_JPsi_Mass2to4p5_v*::HLT_Ele5_SC5_JPsi_Mass2to4p5'),

        # HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*
        cms.InputTag('HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltL1sSingleAndDoubleEGor:HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_L1OR'),
        cms.InputTag('HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter:HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_LeadingLeg'),
        cms.InputTag('HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter:HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_TrailingLeg'),
        cms.InputTag('HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*::HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL'),

        # HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*
        cms.InputTag('HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltL1sSingleAndDoubleEGor:HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_L1OR'),
        cms.InputTag('HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter:HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_LeadingLeg'),
        cms.InputTag('HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter:HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_TrailingLeg'),
        cms.InputTag('HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*::HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL'),

        # HLT_Mu*_TrkIsoVVL_Ele*_CaloIdL_TrackIdL_IsoVL_v*
        cms.InputTag('HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter:HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_ElectronLeg'),
        cms.InputTag('HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter:HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_ElectronLeg'),
        cms.InputTag('HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter:HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_ElectronLeg'),

        # HLT_Mu*_TrkIsoVVL_Ele*_CaloIdL_TrackIdL_IsoVL_DZ_v*
        cms.InputTag('HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter:HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_ElectronLeg'),
        cms.InputTag('HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter:HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_ElectronLeg'),
        cms.InputTag('HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter:HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_ElectronLeg'),
        cms.InputTag('HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter:HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_ElectronLeg'),

        # HLT_*Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v*
        cms.InputTag('HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v*:hltDoubleEle8Mass8Filter:HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_ElectronLeg'),
        cms.InputTag('HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v*:hltElectronMuonInvMassFilter8:HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_ElectronLeg'),

        # single electron triggers
        #cms.InputTag('HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle33CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        #cms.InputTag('HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle18CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle8CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        #cms.InputTag('HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF_v*:hltSingleEle10CaloIdMTrackIdMDphiFilter:HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF_ElectronLeg'),  #Spring15 MC version
        #cms.InputTag('HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_v*:hltSingleEle10CaloIdMTrackIdMDphiFilter:HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_ElectronLeg'),#data version       
        cms.InputTag('HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_v*:hltSingleEle10CaloIdMTrackIdMDphiFilter:HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_ElectronLeg'),
        #cms.InputTag('HLT_Ele33_CaloIdM_TrackIdM_PFJet30_v*:hltEle33CaloIdMGsfTrackIdMDphiFilter:HLT_Ele33_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v*:hltEle23CaloIdMGsfTrackIdMDphiFilter:HLT_Ele23_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),
        #cms.InputTag('HLT_Ele18_CaloIdM_TrackIdM_PFJet30_v*:hltEle18CaloIdMGsfTrackIdMDphiFilter:HLT_Ele18_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*:hltEle17CaloIdMGsfTrackIdMDphiFilter:HLT_Ele17_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v*:hltEle12CaloIdMGsfTrackIdMDphiFilter:HLT_Ele12_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*:hltEle8CaloIdMGsfTrackIdMDphiFilter:HLT_Ele8_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),

        # single electron triggers
        cms.InputTag('HLT_Ele27_eta2p1_WPTight_Gsf_v*::HLT_Ele27_eta2p1_WPTight_Gsf'),
        cms.InputTag('HLT_Ele32_eta2p1_WPTight_Gsf_v*::HLT_Ele32_eta2p1_WPTight_Gsf'),
        cms.InputTag('HLT_Ele105_CaloIdVT_GsfTrkIdT_v*::HLT_Ele105_CaloIdVT_GsfTrkIdT'),
        cms.InputTag('HLT_Ele115_CaloIdVT_GsfTrkIdT_v*::HLT_Ele115_CaloIdVT_GsfTrkIdT') ),
        
        # ### RUN I ###
        # # HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
        # cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltL1sL1DoubleEG137:HLT_Ele17_Ele8_L1sL1DoubleEG137'),
        # cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter:HLT_Ele17_Ele8_LeadingLeg'),
        # cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter:HLT_Ele17_Ele8_TrailingLeg'),
        # cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Ele17_Ele8'),

        # # HLT_Ele27_WP80_v*
        # cms.InputTag('HLT_Ele27_WP80_v*:hltL1sL1SingleEG20ORL1SingleEG22:HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22'),
        # cms.InputTag('HLT_Ele27_WP80_v*::HLT_Ele27_WP80'),

        # # HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*
        # cms.InputTag('HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*:hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter:HLT_Ele17_Ele8_Mass50_LeadingLeg'),
        # cms.InputTag('HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*:hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter:HLT_Ele17_Ele8_Mass50_TrailingLeg'),

        # # HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*
        # cms.InputTag('HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*:hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter:HLT_Ele20_SC4_Mass50_LeadingLeg'),
        # cms.InputTag('HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*:hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter:HLT_Ele20_SC4_Mass50_TrailingLeg'),

        # # HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3
        # cms.InputTag('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*:hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter:HLT_Ele32_SC17_Mass50_LeadingLeg'),
        # cms.InputTag('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*:hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter:HLT_Ele32_SC17_Mass50_TrailingLeg'),

        # # HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
        # # electron related parts only, e.g. electron leg before final filter
        # # and final filter
        # cms.InputTag('HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter:HLT_Mu17_Ele8_TrailingLeg'),
        # cms.InputTag('HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Mu17_Ele8'),

        # # HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
        # # only store final filter for this one, since it is the electron one
        # cms.InputTag('HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Mu8_Ele17')),
        

        processName = cms.untracked.string("HLT"),
        triggerObjectsName = cms.untracked.string("selectedPatTrigger"),
        triggerPrescaleInputTag = cms.untracked.string("patTrigger"),
)
