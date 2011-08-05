import FWCore.ParameterSet.Config as cms

ecalTPFilterMaker = cms.EDProducer('EcalDeadCellEventFilter',

  debug = cms.untracked.bool( False ),

  tpDigiCollection = cms.InputTag("ecalTPSkim"),
  etValToBeFlagged = cms.double(63.75),

  ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
  eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),

  maskedEcalChannelStatusThreshold = cms.int32( 1 ),

  doEEfilter = cms.untracked.bool( True ), # turn it on by default

  makeProfileRoot = cms.untracked.bool( False ),
  profileRootName = cms.untracked.string("deadCellFilterProfile.root" ),

)
