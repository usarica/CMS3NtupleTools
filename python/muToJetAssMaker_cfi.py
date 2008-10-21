import FWCore.ParameterSet.Config as cms

muToJetAssMaker = cms.EDFilter("MuToJetAssMaker",
    # min DR
    minDR = cms.double(0.1)
)


