import FWCore.ParameterSet.Config as cms
from RecoLuminosity.LumiProducer.lumiProducer_cfi import *

luminosityMaker = cms.EDProducer(
  "LuminosityMaker",
  aliasPrefix = cms.untracked.string("ls"),
  lumiSummaryInputTag = cms.InputTag("lumiProducer")
)
