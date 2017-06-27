import FWCore.ParameterSet.Config as cms

muToTrigAssMaker = cms.EDProducer("ObjectToTriggerLegAssMaker",

    aliasPrefix         = cms.untracked.string("mus"),
    objectInputTag      = cms.untracked.InputTag("muonMaker:musp4"),
    cone                = cms.untracked.double(0.2),

    triggers    = cms.untracked.VInputTag(

        ### RUN II ### 
        # Double-lepton matching and single-lepton matching can be done using the TriggerObjects in CMS3, for example,
        # CORE::TriggerSelections::passUnprescaledHLTTrigger(const char* arg, const LorentzVector &obj)               
        # But for LeadingLeg matching we need to know the specific filter passed by the leading leg                   
        # Use:    'TriggerName:FilterName:CMS3name'   --> this requires that muon matches to trigger object passing FilterName of TriggerName 
        # Use:    'TriggerName::CMS3name'             --> this requires that muon matches to trigger object passing the LAST EDFILTER of TriggerName     

        # HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*
        cms.InputTag('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17:HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_LeadingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*:hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8:HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_TrailingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*::HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL'),

        # HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*
        cms.InputTag('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*:hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17:HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_LeadingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*:hltDiMuonGlbFiltered17TrkFiltered8:HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_TrailingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*::HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL'),

        # HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*
        cms.InputTag('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*:hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17:HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_LeadingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*:hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8:HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_TrailingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*::HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ'),

        # HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*
        cms.InputTag('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*:hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17:HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_LeadingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*:hltDiMuonGlbFiltered17TrkFiltered8:HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_TrailingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*::HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ'),

        # HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*
        cms.InputTag('HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*:hltL3fL1sDoubleMu114TkFiltered17Q:HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_LeadingLeg'),
        cms.InputTag('HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*:hltDiTkMuonTkFiltered17TkFiltered8:HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_TrailingLeg'),
        cms.InputTag('HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*::HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ'),

        # HLT_Mu*_TrkIsoVVL_Ele*_CaloIdL_TrackIdL_IsoVL_v*
        cms.InputTag('HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*:hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8:HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_MuonLeg'),
        cms.InputTag('HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered12:HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_MuonLeg'),
        cms.InputTag('HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23:HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_MuonLeg'),

        # HLT_Mu*_TrkIsoVVL_Ele*_CaloIdL_TrackIdL_IsoVL_DZ_v*
        cms.InputTag('HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8:HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_MuonLeg'),
        cms.InputTag('HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered12:HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_MuonLeg'),
        cms.InputTag('HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23:HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_MuonLeg'),
        cms.InputTag('HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23:HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_MuonLeg'),

        # HLT_*Mu8*_Mass8_PFHT300_v*
        cms.InputTag('HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v*:hltElectronMuonInvMassFilter8:HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_MuonLeg'),
        cms.InputTag('HLT_DoubleMu8_Mass8_PFHT300_v*:hltDoubleMu8Mass8L3Filtered:HLT_DoubleMu8_Mass8_PFHT300_MuonLeg'),

        # muon+btag trigger
        cms.InputTag('HLT_Mu10_CentralPFJet30_BTagCSV_p13_v*:hltL3fL1sMu0L1f0L2f3QL3Filtered10Q:HLT_Mu10_CentralPFJet30_BTagCSV_p13_MuonLeg'), 

        # single muon triggers, iso and non-iso, end of 2016
        cms.InputTag('HLT_IsoMu24_v*::HLT_IsoMu24'),
        cms.InputTag('HLT_IsoTkMu24_v*::HLT_IsoTkMu24'),
        cms.InputTag('HLT_Mu50_v*::HLT_Mu50'),
        cms.InputTag('HLT_TkMu50_v*::HLT_TkMu50'),
        cms.InputTag('HLT_Mu55_v*::HLT_Mu55'),

        # fake rate triggers
        cms.InputTag('HLT_Mu8_v*::HLT_Mu8'),
        cms.InputTag('HLT_Mu17_v*::HLT_Mu17'),
        cms.InputTag('HLT_Mu8_TrkIsoVVL_v*::HLT_Mu8_TrkIsoVVL'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_v*::HLT_Mu17_TrkIsoVVL'),

        # single muon triggers 2017 (in addition to some of the above)
        cms.InputTag('HLT_IsoMu27_v*::HLT_IsoMu27'),
    ),
                                  
    processName = cms.untracked.string("HLT"),
    triggerObjectsName = cms.untracked.string("selectedPatTrigger"),
    triggerPrescaleInputTag = cms.untracked.string("patTrigger"),

)
