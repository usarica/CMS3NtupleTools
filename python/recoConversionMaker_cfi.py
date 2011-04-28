import FWCore.ParameterSet.Config as cms

recoConversionMaker = cms.EDProducer("RecoConversionMaker",
                                     recoConversionInputTag = cms.InputTag("allConversions"),
                                     beamSpotInputTag       = cms.InputTag("offlineBeamSpot")
)


