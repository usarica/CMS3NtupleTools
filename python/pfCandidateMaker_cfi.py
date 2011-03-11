import FWCore.ParameterSet.Config as cms

pfCandidateMaker = cms.EDProducer("PFCandidateMaker",
                                  pfCandidatesTag     = cms.InputTag("particleFlow"),
                                  pfElectronsTag      = cms.InputTag("particleFlow","electrons"),
                                  tracksInputTag      = cms.InputTag("generalTracks"),
                                  minDRelectron       = cms.double(0.1)
)
