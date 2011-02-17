from MCProduction2011_cfg import *

process.GlobalTag.globaltag = "START39_V8::All"
process.out = cms.OutputModule(
        "PoolOutputModule",
        dropMetaData = cms.untracked.string("NONE"),
        fileName = cms.untracked.string('ntuple.root')
)
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS2*'))
process.outpath         = cms.EndPath(process.out)
process.p = cms.Path(process.cms2WithEverything)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
