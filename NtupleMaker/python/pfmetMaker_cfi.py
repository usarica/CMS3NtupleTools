import FWCore.ParameterSet.Config as cms


pfmetMaker = cms.EDProducer(
   "PFMETMaker",
   aliasprefix = cms.untracked.string("pfmet"),
   metSrc = cms.InputTag("slimmedMETs"),

   applyMETfix = cms.bool(False),
   metDefaultSrc = cms.InputTag("slimmedMETs"),
   )

pfmetpuppiMaker = cms.EDProducer(
   "PFMETMaker",
   aliasprefix = cms.untracked.string("puppimet"),
   metSrc = cms.InputTag("slimmedMETsPuppi"),

   applyMETfix = cms.bool(False),
   metDefaultSrc = cms.InputTag("slimmedMETsPuppi"),
   )

