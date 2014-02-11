# Import Python Modules
import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff        import *
#from JetMETCorrections.Type1MET.MetType1Corrections_cff import *
from JetMETCorrections.Type1MET.caloMETCorrections_cff import *

# CMS2
process = cms.Process("CMS2")

# Version Control For Python Configuration Files
process.configurationMetadata = cms.untracked.PSet(
        version    = cms.untracked.string('$Revision: 1.11 $'),
        annotation = cms.untracked.string('CMS2'),
        name       = cms.untracked.string('CMS2 test configuration')
)

# load event level configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration/EventContent/EventContent_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("RecoJets.Configuration.RecoJPTJets_cff")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
#process.eleIsoSequence = setupPFElectronIso(process, 'gedGsfElectrons')

process.load("CMS2.NtupleMaker.cms2CoreSequences_cff")
#process.CMS2Reco *= process.pfParticleSelectionSequence
#process.CMS2Reco *= process.eleIsoSequence
process.load("CMS2.NtupleMaker.cms2GENSequence_cff")
process.load('CMS2.NtupleMaker.pixelDigiMaker_cfi')
process.load("CMS2.NtupleMaker.cms2HFCleaningSequence_cff")
process.load("CMS2.NtupleMaker.cms2HcalCleaningSequence_cff")
process.load("CMS2.NtupleMaker.cms2PFSequence_cff")
process.load('RecoJets.Configuration.RecoPFJets_cff') # Import the Jet RECO modules
process.kt6PFJets.doRhoFastjet  = True                # Turn-on the FastJet density calculation
process.ak5PFJets.doAreaFastjet = True                # Turn-on the FastJet jet area calculation for your favorite algorithm


####################
# MET Filters 2012 #
####################

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
#process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')
process.load('RecoMET.METFilters.hcalLaserEventFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellDeltaRFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.load('RecoMET.METFilters.inconsistentMuonPFCandidateFilter_cfi')
process.load('RecoMET.METFilters.greedyMuonPFCandidateFilter_cfi')


#undo what's pulled in by including Reconstruction_cff
#it relies on transient steps introduced in PF in 44X (back-fill)
process.pfPileUp.PFCandidates = cms.InputTag("particleFlowPtrs")
process.pfNoPileUp.bottomCollection = cms.InputTag("particleFlowPtrs") 
process.pfPileUpIso.PFCandidates = cms.InputTag("particleFlowPtrs")
process.pfNoPileUpIso.bottomCollection = cms.InputTag("particleFlowPtrs") 

#
#metJESCorAK5CaloJet.inputUncorJetsLabel = cms.string("ak5CaloJets")


# Input
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(
      #'file:/nfs-3/userdata/cms2/RelValProdTTbar_CMSSW_5_2_2-START52_V4-v2/82900610-FC74-E111-B01A-0018F3D09628.root'   # CMSSW_5_2_2
      #'file:/nfs-3/userdata/cms2/cms2_validation/RelValProdTTbar_CMSSW_5_2_3-START52_V5-v1/144D6226-2C7A-E111-8629-003048678B7C.root'   # CMSSW_5_2_3
      #'file:/nfs-3/userdata/cms2/cms2_validation/RelValTTbar_CMSSW_5_2_3_patch3-START52_V9_special_120410-v1/F8D46BF0-1083-E111-9B6A-001A92811728.root'   # CMSSW_5_2_3_patch3
      #'file:/nfs-3/userdata/fgolf/synchronization/isoSynchFile_DoubleMu191700.root' #muon synch
      #'file:/nfs-3/userdata/fgolf/synchronization/DoubleElectron_Run2012_synchFile.root' #electron synch


      #'file:/nfs-3/userdata/cms2/cms2_validation/SingleMu_CMSSW_5_2_5_cand1-GR_R_52_V7_RelVal_mu2011A-v1_RECO/0CF15A58-7B91-E111-9FDA-002618FDA237.root' # CMSSW_5_2_5
      #'file:/nfs-3/userdata/cms2/cms2_validation/SingleMu_CMSSW_5_3_1-GR_R_53_V2_RelVal_mu2011A-v1_RECO/B4750996-E5A4-E111-8DEC-0030486791F2.root'        # CMSSW_5_3_1

      #'file:../test/TTJets_D6EA05D1-26C2-E111-BF3A-003048FFD7A2.root' # 52x TTJets
        'file:/nfs-3/userdata/jgran/7_0_0_pre12_RelValProdTTbar.root'

   ),
    #--- Uncomment to emulate AOD with RECO --- #
    #inputCommands = process.AODEventContent.outputCommands,
)

# Speed up I/O from castor
process.AdaptorConfig = cms.Service (
  "AdaptorConfig",
  stats = cms.untracked.bool(True),
  enable = cms.untracked.bool(True),
  cacheHint = cms.untracked.string("lazy-download"),
  readHint = cms.untracked.string("auto-detect")
)

# Options
process.options                       = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound'))
#Rethrow = cms.untracked.vstring('ProductNotFound'),
process.source.noEventSort            = cms.untracked.bool( True )
process.MessageLogger.cerr.threshold  = ''

# Number of Events to Process
process.maxEvents                     = cms.untracked.PSet( input = cms.untracked.int32(10) )

# Hypothesis cuts
process.hypDilepMaker.TightLepton_PtCut  = cms.double(20.0)
process.hypDilepMaker.LooseLepton_PtCut  = cms.double(10.0)
process.hypTrilepMaker.TightLepton_PtCut = cms.double(20.0)
process.hypTrilepMaker.LooseLepton_PtCut = cms.double(10.0)

# Event Maker
process.eventMaker.datasetName = cms.string("")
process.eventMaker.CMS2tag     = cms.string("")

