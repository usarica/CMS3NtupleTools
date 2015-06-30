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

        cms.InputTag('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17:HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_LeadingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*::HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL'),

        cms.InputTag('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*:hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17:HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_LeadingLeg'),
        cms.InputTag('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*::HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL'),

        ### RUN I ### 
        # HLT_Mu17_Mu8_v*
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL1sL1DoubleMu10MuOpen:HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen'),
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8:HLT_Mu17_Mu8_TrailingLeg'),
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17:HLT_Mu17_Mu8_LeadingLeg'),
        cms.InputTag('HLT_Mu17_Mu8_v*::HLT_Mu17_Mu8'),

        # HLT_Mu17_TkMu8_v*
        cms.InputTag('HLT_Mu17_TkMu8_v*:hltL2fL1sDoubleMu10MuOpenL1f0L2Filtered10:HLT_Mu17_TkMu8_TrailingLeg'),
        cms.InputTag('HLT_Mu17_TkMu8_v*:hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17:HLT_Mu17_TkMu8_LeadingLeg'),
        cms.InputTag('HLT_Mu17_TkMu8_v*:hltDiMuonGlbFiltered17TrkFiltered8:HLT_Mu17_TkMu8_TrailingLegTrkFiltered'),
        cms.InputTag('HLT_Mu17_TkMu8_v*::HLT_Mu17_TkMu8'),

        # HLT_IsoMu24_eta2p1_v*
        cms.InputTag('HLT_IsoMu24_eta2p1_v*:hltL1sMu16Eta2p1:HLT_IsoMu24_eta2p1_L1sMu16Eta2p1'),
        cms.InputTag('HLT_IsoMu24_eta2p1_v*::HLT_IsoMu24_eta2p1'),

        # HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
        # muon related parts only, e.g. muon leg before final filter
        # and final filter
        cms.InputTag('HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltL1Mu12EG7L3MuFiltered17:HLT_Mu17_Ele8_LeadingLeg'),
        cms.InputTag('HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Mu17_Ele8'),

        # HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
        # as above
        cms.InputTag('HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltL1MuOpenEG12L3Filtered8:HLT_Mu8_Ele17_TrailingLeg'),
        cms.InputTag('HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Mu8_Ele17')),

        processName = cms.untracked.string("HLT"),
        triggerObjectsName = cms.untracked.string("selectedPatTrigger"),

)
