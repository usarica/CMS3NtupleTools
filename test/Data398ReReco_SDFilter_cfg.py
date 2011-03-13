from RecoConfiguration2011_cfg import *

# Global Tag
process.GlobalTag.globaltag = "FT_R_39X_V4A::All"

# Load Filters
process.load("CMS2.NtupleMaker.sdFilter_cfi")
process.filter = cms.Path(process.sdFilter)

# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
        SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('filter') ),
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
process.cms2WithEverything = cms.Sequence( process.sdFilter * process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence )
process.p                  = cms.Path( process.cms2WithEverything )

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(True)
