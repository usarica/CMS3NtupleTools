import FWCore.ParameterSet.Config as cms

eeBadRecovMaker = cms.EDProducer("EEBadRecovMaker",

  EERecHitSource = cms.InputTag("reducedEcalRecHitsEE"),
  MinRecovE = cms.double(30)   
                            
)


