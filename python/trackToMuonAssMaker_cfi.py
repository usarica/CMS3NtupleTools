import FWCore.ParameterSet.Config as cms

trackToMuonAssMaker = cms.EDFilter("TrackToMuonAssMaker",
                                   minDR = cms.double(0.00001),
                                   musInputTag = cms.InputTag("muonMaker"),
                                   trksInputTag = cms.InputTag("trackMaker")
)


