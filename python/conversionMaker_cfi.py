import FWCore.ParameterSet.Config as cms

conversionMaker = cms.EDProducer("ConversionMaker",
    minFracSharedHits = cms.double(0.45)
)


