from CMS2.NtupleMaker.RecoConfiguration2011_44X_cfg import *

# Global Tag
process.GlobalTag.globaltag = "GR_R_44_V4::All"

# Load Filters
# Lower default hypothesis thresholds
# Later we filter dilepton events to have either
# * 20/10 or
# * 5/5 + 100GeV sumJetPt
process.hypDilepMaker.TightLepton_PtCut  = cms.double(5.0)
process.hypDilepMaker.LooseLepton_PtCut  = cms.double(5.0)
process.hypTrilepMaker.TightLepton_PtCut = cms.double(5.0)
process.hypTrilepMaker.LooseLepton_PtCut = cms.double(5.0)

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
process.cms2WithEverything = cms.Sequence( process.ak5PFJets * process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence)
process.p                  = cms.Path( process.cms2WithEverything )

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(True)

process.source.fileNames = [
    'root://xrootd.unl.edu//store/relval/CMSSW_4_4_0_pre9/SingleMu/RECO/GR_R_44_V4_RelVal_mu2011A-v1/0000/EAFCB646-ACD2-E011-BDA0-002618943905.root'
                            ]

process.out.fileName = 'ntuple.root'
process.maxEvents.input = 10
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
#process.Timing =cms.Service("Timing")        

#Slim CMS2
from CMS2.NtupleMaker.SlimCms2_cff import slimcms2
process.out.outputCommands.extend(slimcms2)
