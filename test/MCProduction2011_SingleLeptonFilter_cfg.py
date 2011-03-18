from CMS2.NtupleMaker.RecoConfiguration2011_cfg import *

# Global Tag
process.GlobalTag.globaltag = "MC_311_V2::All"

# Load Filters
process.load('CMS2.NtupleMaker.aSkimFilter_cfi')
process.load('CMS2.NtupleMaker.monolepGenFilter_cfi')

# Single Lepton Filter
process.EventSelectionSingleFilt = cms.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('pWithRecoLepton', 'pWithGenLepton')
  )
)

# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
        process.EventSelectionSingleFilt,
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
process.cms2WithEverything = cms.Sequence( process.ak5PFJets * process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence * process.cms2GENSequence )
process.p                  = cms.Path( process.cms2WithEverything )
process.pWithRecoLepton    = cms.Path(process.cms2WithEverything * process.aSkimFilter   )
process.pWithGenLepton     = cms.Path(process.cms2WithEverything * process.monolepGenFilter  )

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(False)
