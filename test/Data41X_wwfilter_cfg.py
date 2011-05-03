from CMS2.NtupleMaker.RecoConfiguration2011_cfg import *

# Global Tag
process.GlobalTag.globaltag = "GR_R_311_V2::All"

#
process.cms2WithEverything = cms.Sequence( process.ak5PFJets * process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence )

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(True)


# Load Filters
process.load("CMS2.NtupleMaker.wwfilter_cfi")
process.CMS2_Ele10Mu10_IsoIdMET = cms.Path( process.ele10mu10IsoIdMET * process.cms2WithEverything )
process.CMS2_Ele10Mu10_IsoId    = cms.Path( process.ele10mu10IsoId    * process.cms2WithEverything )

# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
        SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('CMS2_Ele10Mu10_IsoIdMET',
                                                                      'CMS2_Ele10Mu10_IsoId') ),
        fileName     = cms.untracked.string('ntuple.root'),
        # dropMetaData = cms.untracked.string("NONE")
)
process.outpath      = cms.EndPath(process.out)

# Branches
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('keep *_ele10*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS2*'))

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring( 'file:/tas/dmytro/tmp/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1_AODSIM_7840115C-5D4F-E011-BD85-001E0B48D9A4.root' )
)
             
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
