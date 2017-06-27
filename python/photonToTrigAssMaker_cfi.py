import FWCore.ParameterSet.Config as cms

photonToTrigAssMaker = cms.EDProducer("ObjectToTriggerLegAssMaker",

    aliasPrefix         = cms.untracked.string("photons"),
    objectInputTag      = cms.untracked.InputTag("photonMaker:photonsp4"),
    cone                = cms.untracked.double(0.2),

    triggers    = cms.untracked.VInputTag(

        ### RUN II ### 
        # Double-lepton matching and single-lepton matching can be done using the TriggerObjects in CMS3, for example,
        # CORE::TriggerSelections::passUnprescaledHLTTrigger(const char* arg, const LorentzVector &obj)               
        # But for LeadingLeg matching we need to know the specific filter passed by the leading leg                   
        # Use:    'TriggerName:FilterName:CMS3name'   --> this requires that photon matches to trigger object passing FilterName of TriggerName 
        # Use:    'TriggerName::CMS3name'             --> this requires that photon matches to trigger object passing the LAST EDFILTER of TriggerName
        
        # single photon triggers
        cms.InputTag('HLT_Photon22_R9Id90_HE10_IsoM_v*::HLT_Photon22_R9Id90_HE10_IsoM'),
        cms.InputTag('HLT_Photon30_R9Id90_HE10_IsoM_v*::HLT_Photon30_R9Id90_HE10_IsoM'),
        cms.InputTag('HLT_Photon36_R9Id90_HE10_IsoM_v*::HLT_Photon36_R9Id90_HE10_IsoM'),
        cms.InputTag('HLT_Photon50_R9Id90_HE10_IsoM_v*::HLT_Photon50_R9Id90_HE10_IsoM'),
        cms.InputTag('HLT_Photon75_R9Id90_HE10_IsoM_v*::HLT_Photon75_R9Id90_HE10_IsoM'),
        cms.InputTag('HLT_Photon90_R9Id90_HE10_IsoM_v*::HLT_Photon90_R9Id90_HE10_IsoM'),
        cms.InputTag('HLT_Photon120_R9Id90_HE10_IsoM_v*::HLT_Photon120_R9Id90_HE10_IsoM'),
        cms.InputTag('HLT_Photon165_R9Id90_HE10_IsoM_v*::HLT_Photon165_R9Id90_HE10_IsoM'),
        cms.InputTag('HLT_Photon165_HE10_v*::HLT_Photon165_HE10'),
        cms.InputTag('HLT_Photon200_v*::HLT_Photon200'),
    ),

    processName = cms.untracked.string("HLT"),
    triggerObjectsName = cms.untracked.string("selectedPatTrigger"),
    triggerPrescaleInputTag = cms.untracked.string("patTrigger"),

)
