import FWCore.ParameterSet.Config as cms

recoConversionMaker = cms.EDProducer("RecoConversionMaker",
                                     recoConversionInputTag = cms.InputTag("reducedEgamma","reducedConversions"),
                                     beamSpotInputTag       = cms.InputTag("offlineBeamSpot")
)


