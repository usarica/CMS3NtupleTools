import FWCore.ParameterSet.Config as cms

pixelDigiMaker = cms.EDFilter("PixelDigiMaker",
    pixelsInputTag          = cms.InputTag("siPixelClusters")
)


