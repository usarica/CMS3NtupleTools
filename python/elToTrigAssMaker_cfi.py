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
        cms.InputTag('HLT_Ele25WP60_Ele8_Mass55_v*:hltEle25WP60Ele8TrackIsoFilter:HLT_Ele25WP60_Ele8_Mass55_LeadingLeg'),
        cms.InputTag('HLT_Ele25WP60_Ele8_Mass55_v*::HLT_Ele25WP60_Ele8_Mass55'),

        cms.InputTag('HLT_Ele25WP60_SC4_Mass55_v*:hltEle25WP60SC4TrackIsoFilter:HLT_Ele25WP60_SC4_Mass55_LeadingLeg'),
        cms.InputTag('HLT_Ele25WP60_SC4_Mass55_v*::HLT_Ele25WP60_SC4_Mass55'),

        cms.InputTag('HLT_Ele5_SC5_JPsi_Mass2to4p5_v*:hltEle5SC5JPsiTrackIsoFilter:HLT_Ele5_SC5_JPsi_Mass2to4p5_LeadingLeg'),
        cms.InputTag('HLT_Ele5_SC5_JPsi_Mass2to4p5_v*::HLT_Ele5_SC5_JPsi_Mass2to4p5'),

        cms.InputTag('HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*:hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter:HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_LeadingLeg'),
        cms.InputTag('HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*::HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL'),

        ### RUN I ###
        # HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltL1sL1DoubleEG137:HLT_Ele17_Ele8_L1sL1DoubleEG137'),
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter:HLT_Ele17_Ele8_LeadingLeg'),
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter:HLT_Ele17_Ele8_TrailingLeg'),
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Ele17_Ele8'),

        # HLT_Ele27_WP80_v*
        cms.InputTag('HLT_Ele27_WP80_v*:hltL1sL1SingleEG20ORL1SingleEG22:HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22'),
        cms.InputTag('HLT_Ele27_WP80_v*::HLT_Ele27_WP80'),

        # HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*
        cms.InputTag('HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*:hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter:HLT_Ele17_Ele8_Mass50_LeadingLeg'),
        cms.InputTag('HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*:hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter:HLT_Ele17_Ele8_Mass50_TrailingLeg'),

        # HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*
        cms.InputTag('HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*:hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter:HLT_Ele20_SC4_Mass50_LeadingLeg'),
        cms.InputTag('HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*:hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter:HLT_Ele20_SC4_Mass50_TrailingLeg'),

        # HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3
        cms.InputTag('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*:hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter:HLT_Ele32_SC17_Mass50_LeadingLeg'),
        cms.InputTag('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*:hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter:HLT_Ele32_SC17_Mass50_TrailingLeg'),

        # HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
        # electron related parts only, e.g. electron leg before final filter
        # and final filter
        cms.InputTag('HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*:hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter:HLT_Mu17_Ele8_TrailingLeg'),
        cms.InputTag('HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Mu17_Ele8'),

        # HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*
        # only store final filter for this one, since it is the electron one
        cms.InputTag('HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Mu8_Ele17')),

        processName = cms.untracked.string("HLT"),
        triggerObjectsName = cms.untracked.string("selectedPatTrigger"),
)
