import FWCore.ParameterSet.Config as cms

jetToMuAssMaker = cms.EDFilter("JetToMuAssMaker",
    # min DR
    minDR = cms.double(0.4)
)


