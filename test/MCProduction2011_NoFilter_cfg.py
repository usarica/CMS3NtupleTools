from RecoConfiguration2011_cfg import *

# Global Tag
process.GlobalTag.globaltag = "START39_V8::All"

# Load Filters

# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
        fileName     = cms.untracked.string('ntuple.root'),
        dropMetaData = cms.untracked.string("NONE")
)
process.outpath      = cms.EndPath(process.out)

# Branches
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS2*'))

#
process.cms2WithEverything = cms.Sequence( process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence * process.cms2GENSequence )
process.p                  = cms.Path( process.cms2WithEverything )

#
process.eventMaker.isData                        = cms.bool(False)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
