from CMS3.NtupleMaker.RecoConfiguration2012_cfg import *

# Global Tag
process.GlobalTag.globaltag = "START53_V7A::All"

# Load Filters
process.load('CMS3.NtupleMaker.aSkimFilter_cfi')
process.load('CMS3.NtupleMaker.monolepGenFilter_cfi')
process.load("CMS3.NtupleMaker.sdFilter_cfi")

process.sdFilter.filterName_=cms.string("Photon")
process.sdFilter.photonJet_dotrig_=cms.bool(False)

#now 3 paths
process.EventSelectionSingleFilt = cms.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('pWithRecoLepton', 'pWithGenLepton', 'pPhoton')
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
process.cms2WithEverything.remove(process.jptMaker)
process.cms2WithEverything.remove(process.hypTrilepMaker)
process.cms2WithEverything.remove(process.hypQuadlepMaker)
process.pWithRecoLepton    = cms.Path(process.cms2WithEverything * process.aSkimFilter   )
process.pWithGenLepton     = cms.Path(process.cms2WithEverything * process.monolepGenFilter  )
process.pPhoton            = cms.Path(process.cms2WithEverything * process.sdFilter)

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(False)
process.luminosityMaker.isData                   = process.eventMaker.isData
