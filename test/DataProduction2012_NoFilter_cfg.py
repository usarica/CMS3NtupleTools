from CMS3.NtupleMaker.RecoConfiguration2015_cfg import *

#Global Tag
process.GlobalTag.globaltag = "GR_R_74_V8"

#Input
process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring('file:/home/users/gzevi/ntupling/CMSSW_7_4_0_pre9/src/CMS3/NtupleMaker/4469D07F-4ADB-E411-B2EB-0025905AA9CC.root')
)

#Output
process.out = cms.OutputModule("PoolOutputModule",
#  fileName     = cms.untracked.string('ntupleDoubleEleData4469D07F.root'),
  fileName     = cms.untracked.string('ntuple.root'),
  dropMetaData = cms.untracked.string("NONE")
)
process.outpath = cms.EndPath(process.out)

#Max Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(51) )

#Branches
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS3*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS3*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS3*'))

#Makers
process.p = cms.Path( 
  process.egmGsfElectronIDSequence *     
  process.beamSpotMaker *
  process.vertexMaker *
  process.secondaryVertexMaker *
  process.pfCandidateMaker *
  process.eventMaker *
  process.electronMaker *
  process.muonMaker *
  process.pfJetMaker *
  process.pfJetPUPPIMaker *
  process.subJetMaker *
  process.pfmetMaker *
  process.hltMakerSequence *
  process.pftauMaker *
  process.photonMaker *
#  process.genMaker *
#  process.genJetMaker *
  process.muToTrigAssMaker *  # requires muonMaker
  process.elToTrigAssMaker *  # requires electronMaker
#  process.candToGenAssMaker * # requires electronMaker, muonMaker, pfJetMaker, photonMaker
#  process.pdfinfoMaker *
  process.puSummaryInfoMaker *
  process.recoConversionMaker *
  process.metFilterMaker *
  process.hcalNoiseSummaryMaker *
  process.miniAODrhoSequence *
  process.hypDilepMaker
)

#Options
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.eventMaker.isData                        = cms.bool(True)
process.luminosityMaker.isData                   = process.eventMaker.isData
process.pfmetMaker.isData                        = process.eventMaker.isData

##Slim CMS3
#from CMS3.NtupleMaker.SlimCms3_cff import slimcms3
#process.out.outputCommands.extend(slimcms3)
