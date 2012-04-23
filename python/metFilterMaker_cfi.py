import FWCore.ParameterSet.Config as cms

metFilterMaker = cms.EDProducer(

  "MetFilterMaker",
  aliasPrefix             = cms.untracked.string("filt"),

  ecalBEInputTag            = cms.InputTag("cms2EcalDeadCellBoundaryEnergyFilter"   ),
  ecalDRInputTag            = cms.InputTag("cms2EcalDeadCellDeltaRFilter"           ),
  ecalTPInputTag            = cms.InputTag("cms2EcalDeadCellTriggerPrimitiveFilter" ),
  greedyMuonInputTag        = cms.InputTag("cms2greedyMuonPFCandidateFilter"        ),
  hcalLaserEventInputTag    = cms.InputTag("cms2hcalLaserEventFilter"               ),
  inconsistentMuonInputTag  = cms.InputTag("cms2inconsistentMuonPFCandidateFilter"  ),
  jetIDFailureInputTag      = cms.InputTag("cms2jetIDFailureFilter"                 ),
  multiEventFailureInputTag = cms.InputTag("cms2multiEventFailureFilter"            ),
  trackingFailureInputTag   = cms.InputTag("cms2trackingFailureFilter"              )

  #ecalBEInputTag            = cms.InputTag("cms2EcalDeadCellBoundaryEnergyFilter"   , "EcalDeadCellBoundaryEnergyFilter"  ),
  #ecalDRInputTag            = cms.InputTag("cms2EcalDeadCellDeltaRFilter"           , "EcalDeadCellDeltaRFilter"          ),
  #ecalTPInputTag            = cms.InputTag("cms2EcalDeadCellTriggerPrimitiveFilter" , "EcalDeadCellTriggerPrimitiveFilter"),
  #greedyMuonInputTag        = cms.InputTag("cms2greedyMuonPFCandidateFilter"        , "greedyMuonPFCandidateFilter"       ),
  #hcalLaserEventInputTag    = cms.InputTag("cms2hcalLaserEventFilter"               , "hcalLaserEventFilter"              ),
  #inconsistentMuonInputTag  = cms.InputTag("cms2inconsistentMuonPFCandidateFilter"  , "inconsistentMuonPFCandidateFilter" ),
  #jetIDFailureInputTag      = cms.InputTag("cms2jetIDFailureFilter"                 , "jetIDFailureFilter"                ),
  #multiEventFailureInputTag = cms.InputTag("cms2multiEventFailureFilter"            , "multiEventFailureFilter"           ),
  #trackingFailureInputTag   = cms.InputTag("cms2trackingFailureFilter"              , "trackingFailureFilter"             )

)
