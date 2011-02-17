from MCProduction2011_cfg import *

# Single Lepton Filter
process.GlobalTag.globaltag = "START39_V8::All"
process.EventSelectionSingleFilt = cms.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('pWithRecoLepton', 'pWithGenLepton')
  )
)
process.out = cms.OutputModule(
        "PoolOutputModule",
        process.EventSelectionSingleFilt,
        dropMetaData = cms.untracked.string("NONE"),
        fileName = cms.untracked.string('ntuple.root')
)
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS2*'))
process.outpath         = cms.EndPath(process.out)
process.load('CMS2.NtupleMaker.aSkimFilter_cfi')
process.load('CMS2.NtupleMaker.monolepGenFilter_cfi')
process.p               = cms.Path(process.cms2WithEverything)
process.pWithRecoLepton = cms.Path(process.cms2WithEverything * process.aSkimFilter   )
process.pWithGenLepton  = cms.Path(process.cms2WithEverything * process.monolepGenFilter  )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
