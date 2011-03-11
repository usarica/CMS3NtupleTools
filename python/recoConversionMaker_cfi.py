import FWCore.ParameterSet.Config as cms

recoConversionMaker = cms.EDProducer("RecoConversionMaker",
                            recoConversionInputTag = cms.InputTag("trackerOnlyConversions"),
                            beamSpotInputTag       = cms.InputTag("offlineBeamSpot")
)


