import FWCore.ParameterSet.Config as cms

pfCandidateMaker = cms.EDProducer("PFCandidateMaker",
                                  pfCandidatesTag     = cms.InputTag("particleFlow"),
                                  tracksInputTag      = cms.InputTag("generalTracks")
)
