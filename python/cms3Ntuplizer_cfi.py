import FWCore.ParameterSet.Config as cms

cms3ntuple = cms.EDAnalyzer(
   "CMS3Ntuplizer",

   year = cms.int32(-1), # Must be overriden by main_pset
   treename = cms.untracked.string("Events"),

   isMC = cms.bool(False),

   electronSrc   = cms.InputTag("electronMaker"),
   photonSrc   = cms.InputTag("photonMaker"),
   muonSrc   = cms.InputTag("muonMaker"),

   genInfoSrc = cms.InputTag("genMaker"),

)
