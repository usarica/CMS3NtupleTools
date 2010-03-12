import FWCore.ParameterSet.Config as cms

conversionMaker = cms.EDProducer("ConversionMaker",
	aliasPrefix = cms.untracked.string("trks"),
    minFracSharedHits = cms.double(0.45)
)


