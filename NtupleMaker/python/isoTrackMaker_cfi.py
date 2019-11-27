import FWCore.ParameterSet.Config as cms

isoTrackMaker = cms.EDProducer(
   "IsoTrackMaker",

   aliasprefix      = cms.untracked.string("isotracks"),
   year = cms.int32(-1), # Must be overriden by main_pset

   isoTracksTag        = cms.InputTag("isolatedTracks"),
   lostTracksTag       = cms.InputTag("lostTracks"),
   pfCandidatesTag     = cms.InputTag("packedPFCandidates"),
   )
