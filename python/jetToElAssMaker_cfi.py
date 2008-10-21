import FWCore.ParameterSet.Config as cms

jetToElAssMaker = cms.EDFilter("JetToElAssMaker",
    # min DR
    minDR = cms.double(0.1)
)


