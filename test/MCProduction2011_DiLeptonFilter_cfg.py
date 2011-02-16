from MCProduction2011_cfg import *

# Dilepton Filter
process.hypDilepMaker.TightLepton_PtCut  = cms.double(20.0)
process.hypDilepMaker.LooseLepton_PtCut  = cms.double(10.0)
process.hypTrilepMaker.TightLepton_PtCut = cms.double(20.0)
process.hypTrilepMaker.LooseLepton_PtCut = cms.double(10.0)
process.EventSelectionDilFilt = cms.PSet (
  SelectEvents = cms.untracked.PSet (
    SelectEvents = cms.vstring('pWithHyp', 'pWithGenHyp')
  )
)
process.out = cms.OutputModule(
        "PoolOutputModule",
        process.EventSelectionDilFilt,
        dropMetaData = cms.untracked.string("NONE"),
        fileName = cms.untracked.string('ntuple.root')
)
process.load("CMS2.NtupleMaker.hypFilter_cfi")
process.load("CMS2.NtupleMaker.dilepGenFilter_cfi")
process.pWithHyp    = cms.Path(process.cms2WithEverything * process.hypFilter)
process.pWithGenHyp = cms.Path(process.cms2WithEverything * process.dilepGenFilter)
