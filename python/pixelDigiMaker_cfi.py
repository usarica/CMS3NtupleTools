import FWCore.ParameterSet.Config as cms

pixelDigiMaker = cms.EDProducer("PixelDigiMaker",
	aliasPrefix = cms.untracked.string("pxl"),
    pixelsInputTag          = cms.InputTag("siPixelClusters")
)


