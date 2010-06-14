import FWCore.ParameterSet.Config as cms

trackToMuonAssMaker = cms.EDProducer("TrackToMuonAssMaker",
	aliasPrefix = cms.untracked.string("trk"),
                                   minDR = cms.double(0.00001),
                                   musInputTag = cms.InputTag("muonMaker"),
                                   trksInputTag = cms.InputTag("trackMaker")
)


