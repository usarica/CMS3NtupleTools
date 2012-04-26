import FWCore.ParameterSet.Config as cms

muToTrigAssMaker = cms.EDProducer("ObjectToTriggerLegAssMaker",

    aliasPrefix         = cms.untracked.string("mus"),
    objectInputTag      = cms.untracked.InputTag("muonMaker:musp4"),
    cone                = cms.untracked.double(0.2),

    triggers    = cms.untracked.VInputTag(

        # HLT_Mu17_Mu8_v*
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL1sL1DoubleMu10MuOpen:HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen'),
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8:HLT_Mu17_Mu8_TrailingLeg'),
        cms.InputTag('HLT_Mu17_Mu8_v*:hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17:HLT_Mu17_Mu8_LeadingLeg'),
        cms.InputTag('HLT_Mu17_Mu8_v*::HLT_Mu17_Mu8'))

)

