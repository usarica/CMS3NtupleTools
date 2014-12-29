from CMS3.NtupleMaker.RecoConfiguration2015_cfg import *

#Global Tag
process.GlobalTag.globaltag = "PHYS14_25_V1::All"

#Input
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring('file:/home/users/namin/stop/cms3/CMSSW_7_2_0/src/SMS-T1tttt_PU20bx25_tsg_PHYS14_25_V1.root')
)

#Output
process.out = cms.OutputModule("PoolOutputModule",
  fileName     = cms.untracked.string('ntuple.root'),
  dropMetaData = cms.untracked.string("NONE")
)
process.outpath = cms.EndPath(process.out)

#Max Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(104) )

#Branches
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS2*'))

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

#Slim CMS2
from CMS3.NtupleMaker.SlimCms2_cff import slimcms2
process.out.outputCommands.extend(slimcms2)
