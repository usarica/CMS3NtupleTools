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
        cms.InputTag('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*::HLT_Ele17_Ele8'))

)

