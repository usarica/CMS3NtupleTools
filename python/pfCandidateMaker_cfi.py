import FWCore.ParameterSet.Config as cms

pfCandidateMaker = cms.EDProducer("PFCandidateMaker",
                                  #pfCandidatesTag     = cms.InputTag("particleFlow"),
                                  pfCandidatesTag     = cms.InputTag("packedPFCandidates"),
                                  #pfElectronsTag      = cms.InputTag("particleFlow","electrons"),
                                  #tracksInputTag      = cms.InputTag("generalTracks"),
                                  #vertexInputTag      = cms.InputTag("offlinePrimaryVertices"),
                                  #minDRelectron       = cms.double(0.1)
)
