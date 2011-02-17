from MCProduction2011_cfg import *

# Dilepton Filter
process.GlobalTag.globaltag = "START39_V8::All"
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
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS2*'))
process.outpath         = cms.EndPath(process.out)
process.load("CMS2.NtupleMaker.hypFilter_cfi")
process.load("CMS2.NtupleMaker.dilepGenFilter_cfi")
process.p           = cms.Path(process.cms2WithEverything)
process.pWithHyp    = cms.Path(process.cms2WithEverything * process.hypFilter)
process.pWithGenHyp = cms.Path(process.cms2WithEverything * process.dilepGenFilter)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
