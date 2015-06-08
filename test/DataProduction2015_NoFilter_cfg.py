import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff        import *

#CMS3
process = cms.Process("CMS3")

#Version Control For Python Configuration Files
process.configurationMetadata = cms.untracked.PSet(
        version    = cms.untracked.string('$Revision: 1.11 $'),
        annotation = cms.untracked.string('CMS3'),
        name       = cms.untracked.string('CMS3 test configuration')
)

# load event level configurations
process.load('Configuration/EventContent/EventContent_cff')
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

#services
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.GlobalTag.globaltag = "GR_R_74_V12"
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold  = ''
process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True),SkipEvent = cms.untracked.vstring('ProductNotFound') )

#Output
process.out = cms.OutputModule("PoolOutputModule",
  fileName     = cms.untracked.string('ntuple.root'),
  dropMetaData = cms.untracked.string("NONE")
)
process.outpath = cms.EndPath(process.out)

#Branches 
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS3*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS3*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS3*'))

#load cff and third party tools
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from JMEAnalysis.JetToolbox.jetToolbox_cff import *
#from JetMETCorrections.Configuration.JetCorrectionProducers_cff import *
from JetMETCorrections.Configuration.CorrectedJetProducersDefault_cff import *
from JetMETCorrections.Configuration.CorrectedJetProducers_cff import *
from JetMETCorrections.Configuration.CorrectedJetProducersAllAlgos_cff import *
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from RecoJets.JetProducers.fixedGridRhoProducerFastjet_cfi import *
process.fixedGridRhoFastjetAll = fixedGridRhoFastjetAll.clone(pfCandidatesTag = 'packedPFCandidates')

#Electron Identification for PHYS 14#
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons',"","PAT")
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#Load Ntuple producer cff
process.load("CMS3.NtupleMaker.cms3CoreSequences_cff")
process.load("CMS3.NtupleMaker.cms3PFSequence_cff")

# Hypothesis cuts
process.hypDilepMaker.TightLepton_PtCut  = cms.double(10.0)
process.hypDilepMaker.LooseLepton_PtCut  = cms.double(10.0)

Options for Input
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring('file:/nfs-7/userdata/jgran/74x_sync/1294BDDB-B7FE-E411-8028-002590596490.root')
)
process.source.noEventSort = cms.untracked.bool( True )

#Max Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(51) )

#Run jet tool box
jetToolbox( process, 'ak4', 'ak4JetSubs', 'out',PUMethod='',miniAOD=True,JETCorrLevels=['L1FastJet','L2Relative', 'L3Absolute'])

process.p = cms.Path( 
  process.metFilterMaker *
  process.hcalNoiseSummaryMaker *
  process.egmGsfElectronIDSequence *     
  process.beamSpotMaker *
  process.vertexMaker *
  process.secondaryVertexMaker *
  process.eventMaker *
  process.pfCandidateMaker *
  process.isoTrackMaker *
  process.recoConversionMaker *
  process.electronMaker *
  process.muonMaker *
  process.pfJetMaker *
  process.pfJetPUPPIMaker*
  process.ak4JetMaker *
  process.subJetMaker *
#  process.ca12subJetMaker *
  process.pfmetMaker *
  process.hltMakerSequence *
  process.pftauMaker *
  process.photonMaker *
  #process.genMaker *
  #process.genJetMaker *
  process.muToTrigAssMaker *  # requires muonMaker
  process.elToTrigAssMaker *  # requires electronMaker
  #process.candToGenAssMaker * # requires electronMaker, muonMaker, pfJetMaker, photonMaker
  #process.pdfinfoMaker *
  #process.puSummaryInfoMaker *
  process.recoConversionMaker *
  process.miniAODrhoSequence *
  process.hypDilepMaker
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.eventMaker.isData                        = cms.bool(True)
process.pfmetMaker.isData                        = process.eventMaker.isData
