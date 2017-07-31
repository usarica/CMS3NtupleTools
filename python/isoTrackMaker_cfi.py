import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

isoTrackMaker = cms.EDProducer("IsoTrackMaker",
                                  pfCandidatesTag     = cms.InputTag("packedPFCandidates","",configProcessName.name),
                                  isoTracksTag        = cms.InputTag("isolatedTracks","",configProcessName.name),
                                  pT_cut       = cms.double(5.0),   ## miniAOD has cuts of 5.0/20.0, but can make them stricter here
                                  pT_cut_noIso = cms.double(20.0)
)
