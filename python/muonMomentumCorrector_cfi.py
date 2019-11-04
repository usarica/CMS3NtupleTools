import FWCore.ParameterSet.Config as cms


correctedMuons = cms.EDProducer(
   "RochesterPATMuonCorrector",

   identifier = cms.string("RoccoR2017"),
   isMC = cms.bool(False),
   src = cms.InputTag("slimmedMuons"),

)

