from CMS2.NtupleMaker.RecoConfiguration2012_cfg import *

# Global Tag
process.GlobalTag.globaltag = "START52_V4A::All"

# Load Filters
process.load('CMS2.NtupleMaker.aSkimFilter_cfi')        # from single lepton filter
process.load('CMS2.NtupleMaker.monolepGenFilter_cfi')   # from single lepton filter
process.load("CMS2.NtupleMaker.hypFilter_cfi")          # from di lepton filter
process.load("CMS2.NtupleMaker.dilepGenFilter_cfi")     # from di lepton filter

# Dilepton Filter
process.EventSelectionSingleOrDilFilt = cms.PSet (    # what does changing this do?
  SelectEvents = cms.untracked.PSet (
    SelectEvents = cms.vstring('pWithRecoLepton', 'pWithGenLepton',        # from single lepton filter
                               'pDiLepton', 'pMultiLepton', 'pWithGenHyp'  #from di lepton filter
                               ) #pWithGenLepton and pWithGenHyp should be redundant, leaving them both in to be safe
  )
)

# Lower default hypothesis thresholds
process.hypDilepMaker.TightLepton_PtCut  = cms.double(5.0)  # from di lepton filter
process.hypDilepMaker.LooseLepton_PtCut  = cms.double(5.0)  # from di lepton filter
process.hypTrilepMaker.TightLepton_PtCut = cms.double(5.0)  # from di lepton filter
process.hypTrilepMaker.LooseLepton_PtCut = cms.double(5.0)  # from di lepton filter

# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
        process.EventSelectionSingleOrDilFilt,
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
process.pWithRecoLepton    = cms.Path( process.cms2WithEverything * process.aSkimFilter   )       # from single lepton filter
process.pWithGenLepton     = cms.Path( process.cms2WithEverything * process.monolepGenFilter  )   # from single lepton filter
process.pDiLepton          = cms.Path( process.cms2WithEverything * process.hypDiLeptonFilter )   # from di lepton filter
process.pMultiLepton       = cms.Path( process.cms2WithEverything * process.hypOtherFilter )      # from di lepton filter
process.pWithGenHyp        = cms.Path( process.cms2WithEverything * process.dilepGenFilter )      # from di lepton filter

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(False)
process.luminosityMaker.isData                   = process.eventMaker.isData
#process.source.fileNames = [
#    #'root://xrootd.unl.edu//store/mc/Summer12/TTJets_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S7_START52_V5-v1/0000/920C7062-4D81-E111-A036-001A92810AE4.root'
#    'file:/nfs-3/userdata/cms2/cms2_validation/RelValTTbar_CMSSW_5_2_3_patch3-START52_V9_special_120410-v1/F8D46BF0-1083-E111-9B6A-001A92811728.root'
#    ]
#process.maxEvents.input = 10
