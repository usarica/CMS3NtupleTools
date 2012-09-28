from CMS2.NtupleMaker.RecoConfiguration2012_cfg import *

# Global Tag
process.GlobalTag.globaltag = "GR_R_52_V7::All"

#
process.cms2WithEverything = cms.Sequence( process.ak5PFJets * process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence )

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(True)
process.luminosityMaker.isData                   = process.eventMaker.isData

# Load Filters
process.load("CMS2.NtupleMaker.wwfilter_cfi")
process.CMS2_Ele10Mu10_NoIsoIdMET  = cms.Path( process.ele10mu10NoIsoIdMET * process.cms2WithEverything )
process.CMS2_Ele10Mu10_NoIsoId     = cms.Path( process.ele10mu10NoIsoId    * process.cms2WithEverything )
process.CMS2_Ele10Mu10_NoIsoIdPt45 = cms.Path( process.ele10mu10NoIsoIdPt45 * process.cms2WithEverything )
process.load("CMS2.NtupleMaker.zzfilter_cfi")
process.CMS2_4L = cms.Path( process.fourLeptons    * process.cms2WithEverything )

# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
        SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('CMS2_Ele10Mu10_NoIsoIdMET',
                                                                      'CMS2_Ele10Mu10_NoIsoId',
                                                                      'CMS2_Ele10Mu10_NoIsoIdPt45',
                                                                      'CMS2_4L') ),
        fileName     = cms.untracked.string('ntuple.root'),
)
process.outpath      = cms.EndPath(process.out)

# Branches
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('keep *_ele10*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('keep *_fourLeptons_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS2*'))

#process.source.fileNames = [
#    'file:/smurf/cerati/Run2012A_DoubleElectron_AOD_PromptReco-v1_000_190_733_86B0E544-E283-E111-9A36-BCAEC53296F4.root',
#    'file:/smurf/cerati/Run2012A_DoubleElectron_AOD_PromptReco-v1_000_191_247_04825687-3588-E111-82CE-BCAEC518FF63.root'
#                            ]
             
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#Slim CMS2
from CMS2.NtupleMaker.SlimCms2_cff import slimcms2
process.out.outputCommands.extend(slimcms2)
