import FWCore.ParameterSet.Config as cms

elToTrigAssMaker = cms.EDProducer("ObjectToTriggerLegAssMaker",

    aliasPrefix         = cms.untracked.string("els"),
    objectInputTag      = cms.untracked.InputTag("electronMaker:elsp4"),
    cone                = cms.untracked.double(0.2),

    triggers    = cms.untracked.VInputTag(

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
        cms.InputTag('HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Mu8_Ele17'))

)
