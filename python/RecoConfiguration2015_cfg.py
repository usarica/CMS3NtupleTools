# Import Python Modules
import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff        import *

# CMS3
process = cms.Process("CMS3")

# Version Control For Python Configuration Files
process.configurationMetadata = cms.untracked.PSet(
        version    = cms.untracked.string('$Revision: 1.11 $'),
        annotation = cms.untracked.string('CMS3'),
        name       = cms.untracked.string('CMS3 test configuration')
)

# load event level configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration/EventContent/EventContent_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("RecoJets.Configuration.RecoJPTJets_cff")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

#Identification for PHYS 14
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("CMS3.NtupleMaker.cms2CoreSequences_cff")
process.load("CMS3.NtupleMaker.cms2GENSequence_cff")
#process.load('CMS3.NtupleMaker.pixelDigiMaker_cfi')
process.load("CMS3.NtupleMaker.cms2PFSequence_cff")
process.load('RecoJets.Configuration.RecoPFJets_cff') # Import the Jet RECO modules
process.kt6PFJets.doRhoFastjet  = True                # Turn-on the FastJet density calculation
process.ak5PFJets.doAreaFastjet = True                # Turn-on the FastJet jet area calculation for your favorite algorithm


####################
# MET Filters 2012 #
####################

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.load('RecoMET.METFilters.hcalLaserEventFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellDeltaRFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.load('RecoMET.METFilters.inconsistentMuonPFCandidateFilter_cfi')
process.load('RecoMET.METFilters.greedyMuonPFCandidateFilter_cfi')


#undo what's pulled in by including Reconstruction_cff
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
        'file:/nfs-3/userdata/jgran/7_0_0_pre12_RelValProdTTbar.root'

   ),
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
process.maxEvents                     = cms.untracked.PSet( input = cms.untracked.int32(50) )

# Hypothesis cuts
process.hypDilepMaker.TightLepton_PtCut  = cms.double(10.0)
process.hypDilepMaker.LooseLepton_PtCut  = cms.double(10.0)

# Event Maker
process.eventMaker.datasetName = cms.string("")
process.eventMaker.CMS3tag     = cms.string("")

