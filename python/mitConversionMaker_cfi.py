#import FWCore.ParameterSet.Config as cms
#
#mitConversionMaker = cms.EDFilter("MITConversionMaker",
#                                  elsInputTag = cms.InputTag("gsfElectrons"),
#                                  mitConversionsInputTag = cms.InputTag("mvfConversionRemoval"),
#                                  ctfTrksInputTag   = cms.InputTag("generalTracks"),
#                                  beamSpotInputTag       = cms.InputTag("offlineBeamSpot")
#)
