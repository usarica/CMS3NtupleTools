import FWCore.ParameterSet.Config as cms

trackToMuonAssMaker = cms.EDFilter("TrackToMuonAssMaker",
    # min DR
    minDR = cms.double(0.00001)
)


