from CMS3.NtupleMaker.RecoConfiguration2015_cfg import *

#Global Tag
process.GlobalTag.globaltag = "MCRUN2_74_V9A"

#Input
process.source = cms.Source("PoolSource",
#fileNames = cms.untracked.vstring('file:///home/users/fgolf/devel/CMSSW_7_4_1/src/CMS3/NtupleMaker/test/AA2CC414-DCEC-E411-B459-0025905A608A.root')
fileNames = cms.untracked.vstring('file:///home/users/gzevi/ntupling/CMSSW_7_4_1/src/CMS3/NtupleMaker/QCD_Pt_1400to1800_Asympt50ns_MCRUN2_74_V9A-v1.root')
)

#Output
process.out = cms.OutputModule("PoolOutputModule",
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
  process.genMaker *
  process.genJetMaker *
  process.muToTrigAssMaker *  # requires muonMaker
  process.elToTrigAssMaker *  # requires electronMaker
  process.candToGenAssMaker * # requires electronMaker, muonMaker, pfJetMaker, photonMaker
  process.pdfinfoMaker *
  process.puSummaryInfoMaker *
  process.recoConversionMaker *
  process.metFilterMaker *
  process.hcalNoiseSummaryMaker *
  process.miniAODrhoSequence *
  process.hypDilepMaker
)

#Options
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.eventMaker.isData                        = cms.bool(False)
process.luminosityMaker.isData                   = process.eventMaker.isData
process.pfmetMaker.isData                        = process.eventMaker.isData

##Slim CMS3
#from CMS3.NtupleMaker.SlimCms3_cff import slimcms3
#process.out.outputCommands.extend(slimcms3)
