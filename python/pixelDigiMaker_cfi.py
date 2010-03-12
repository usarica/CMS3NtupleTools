import FWCore.ParameterSet.Config as cms

pixelDigiMaker = cms.EDFilter("PixelDigiMaker",
	aliasPrefix = cms.untracked.string("pxl"),
    pixelsInputTag          = cms.InputTag("siPixelClusters")
)


