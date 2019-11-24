import FWCore.ParameterSet.Config as cms

lumiFilter = cms.EDFilter("LumiFilter",
  lumisToProcess = cms.untracked(cms.VLuminosityBlockRange()),
)


