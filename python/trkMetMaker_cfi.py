import FWCore.ParameterSet.Config as cms

trkMetMaker = cms.EDProducer("TrkMETMaker",
                             aliasPrefix = cms.untracked.string("trk"),
                             pfcandInputTag     = cms.InputTag("pfCandidateMaker"),
                             hypInputTag        = cms.InputTag("hypDilepMaker"),
                             vertexInputTag     = cms.InputTag("vertexMaker"),
                             dzcut              = cms.double(0.1),
                             drcut              = cms.double(0.1),
                             correctJet         = cms.bool(False)
                             )
