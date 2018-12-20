import FWCore.ParameterSet.Config as cms

isoTrackMaker = cms.EDProducer("IsoTrackMaker",
                                  pfCandidatesTag     = cms.InputTag("packedPFCandidates"),
                                  lostTracksTag       = cms.InputTag("lostTracks"),
                                  isoTracksTag        = cms.InputTag("isolatedTracks"),
                                  pT_cut       = cms.double(5.0),   ## miniAOD has cuts of 5.0/20.0, but can make them stricter here
                                  pT_cut_noIso = cms.double(20.0)
)
