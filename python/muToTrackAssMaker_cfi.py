import FWCore.ParameterSet.Config as cms

muToTrackAssMaker = cms.EDFilter("MuToTrackAssMaker",
    # min DR
    minDR = cms.double(0.00001)
)


