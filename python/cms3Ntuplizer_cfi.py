import FWCore.ParameterSet.Config as cms

cms3ntuplizer = cms.EDAnalyzer(
   "CMS3Ntuplizer",
   treename = cms.untracked.string("Events")
)
