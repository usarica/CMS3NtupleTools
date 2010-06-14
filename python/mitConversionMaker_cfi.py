import FWCore.ParameterSet.Config as cms

mitConversionMaker = cms.EDFilter("MITConversionMaker",
                                  elsInputTag = cms.InputTag("gsfElectrons"),
                                  mitConversionsTag = cms.InputTag("mvfConversionRemoval"),
                                  ctfTrksInputTag   = cms.InputTag("generalTracks")
)
