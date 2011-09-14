from CMS2.NtupleMaker.RecoConfiguration2011_44X_cfg import *

# Global Tag
process.GlobalTag.globaltag = "MC_44_V4::All"

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

# Lower default hypothesis thresholds
# Later we filter dilepton events to have either
# * 20/10 or
# * 5/5 + 100GeV sumJetPt
process.hypDilepMaker.TightLepton_PtCut  = cms.double(5.0)
process.hypDilepMaker.LooseLepton_PtCut  = cms.double(5.0)
process.hypTrilepMaker.TightLepton_PtCut = cms.double(5.0)
process.hypTrilepMaker.LooseLepton_PtCut = cms.double(5.0)

#
process.cms2WithEverything = cms.Sequence( process.ak5PFJets * process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence * process.cms2GENSequence )
process.p                  = cms.Path( process.cms2WithEverything )

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(False)

#process.source.fileNames = [
#    'root://xrootd.unl.edu//store/relval/CMSSW_4_4_0_pre9/RelValSingleMuPt100/GEN-SIM-RECO/START44_V4-v1/0174/ECEF9E0F-A6D3-E011-A44B-00261894380A.root',
#    'root://xrootd.unl.edu//store/relval/CMSSW_4_4_0_pre9/RelValSingleMuPt100/GEN-SIM-RECO/START44_V4-v1/0172/000FD6DC-97D2-E011-9EA7-00261894390C.root'
#                            ]
#process.out.fileName = 'RelValSingleMuPt100_CMSSW_4_4_0_pre9_ntuple.root'
#process.maxEvents.input = -1
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
#process.Timing =cms.Service("Timing")        
