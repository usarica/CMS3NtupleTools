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
        cms.InputTag('HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v*:hltDoubleEle8CaloIdMGsfTrackIdMDphiFilter:HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_ElectronLeg'), # note: actual last electron filter, hltDoubleEle8Mass8Filter, appears to not have info saved for some reason
        cms.InputTag('HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v*:hltElectronMuonInvMassFilter8:HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_ElectronLeg'),

        # single electron fake rate triggers
        cms.InputTag('HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*:hltEle8CaloIdLTrackIdLIsoVLTrackIsoFilter:HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_v*:hltSingleEle10CaloIdMTrackIdMDphiFilter:HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_ElectronLeg'),
        cms.InputTag('HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v*:hltEle23CaloIdMGsfTrackIdMDphiFilter:HLT_Ele23_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*:hltEle17CaloIdMGsfTrackIdMDphiFilter:HLT_Ele17_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v*:hltEle12CaloIdMGsfTrackIdMDphiFilter:HLT_Ele12_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),
        cms.InputTag('HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*:hltEle8CaloIdMGsfTrackIdMDphiFilter:HLT_Ele8_CaloIdM_TrackIdM_PFJet30_ElectronLeg'),

        # single electron triggers
        cms.InputTag('HLT_Ele27_eta2p1_WPTight_Gsf_v*::HLT_Ele27_eta2p1_WPTight_Gsf'),
        cms.InputTag('HLT_Ele32_eta2p1_WPTight_Gsf_v*::HLT_Ele32_eta2p1_WPTight_Gsf'),
        cms.InputTag('HLT_Ele105_CaloIdVT_GsfTrkIdT_v*::HLT_Ele105_CaloIdVT_GsfTrkIdT'),
        cms.InputTag('HLT_Ele115_CaloIdVT_GsfTrkIdT_v*::HLT_Ele115_CaloIdVT_GsfTrkIdT'),

        # single electron triggers 2017
        cms.InputTag('HLT_Ele27_WPTight_Gsf_v*::HLT_Ele27_WPTight_Gsf'),
        cms.InputTag('HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*::HLT_Ele32_WPTight_Gsf_L1DoubleEG'),
        cms.InputTag('HLT_Ele35_WPTight_Gsf_v*::HLT_Ele35_WPTight_Gsf'),
        cms.InputTag('HLT_Ele38_WPTight_Gsf_v*::HLT_Ele38_WPTight_Gsf'),
        cms.InputTag('HLT_Ele40_WPTight_Gsf_v*::HLT_Ele40_WPTight_Gsf'),
    ),
                                  
    processName = cms.untracked.string("HLT"),
    triggerObjectsName = cms.untracked.string("selectedPatTrigger"),
    triggerPrescaleInputTag = cms.untracked.string("patTrigger"),
)
