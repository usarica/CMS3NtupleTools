import FWCore.ParameterSet.Config as cms

metFilterMaker = cms.EDProducer(
  "MetFilterMaker",
  aliasPrefix             = cms.untracked.string("filt"),

  #ecalTPInputTag          = cms.InputTag("ecalDeadCellTriggerPrimitiveFilter"),
  #trackingFailureInputTag = cms.InputTag("trackingFailureFilter")
  ecalTPInputTag          = cms.InputTag("cms2EcalDeadCellTriggerPrimitiveFilter", "ecalDeadCellTriggerPrimitiveFilter"),
  trackingFailureInputTag = cms.InputTag("cms2trackingFailureFilter"             , "trackingFailureFilter"             )


)
