# Import Python Modules
import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff        import *
from JetMETCorrections.Type1MET.MetType1Corrections_cff import *

# CMS2
process = cms.Process("CMS2")

# Version Control For Python Configuration Files
process.configurationMetadata = cms.untracked.PSet(
        version    = cms.untracked.string('$Revision: 1.1 $'),
        annotation = cms.untracked.string('CMS2'),
        name       = cms.untracked.string('CMS2 test configuration')
)

# load event level configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration/EventContent/EventContent_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("RecoJets.Configuration.RecoJPTJets_cff")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")
process.load("CMS2.NtupleMaker.cms2CoreSequences_cff")
process.load("CMS2.NtupleMaker.cms2GENSequence_cff")
process.load('CMS2.NtupleMaker.pixelDigiMaker_cfi')
process.load("CMS2.NtupleMaker.cms2HFCleaningSequence_cff")
process.load("CMS2.NtupleMaker.cms2HcalCleaningSequence_cff")
process.load("CMS2.NtupleMaker.cms2PFSequence_cff")
#
process.load('RecoJets.Configuration.RecoPFJets_cff') # Import the Jet RECO modules
process.kt6PFJets.doRhoFastjet  = True                # Turn-on the FastJet density calculation
process.ak5PFJets.doAreaFastjet = True                # Turn-on the FastJet jet area calculation for your favorite algorithm
#
#process.ak5CaloL1Offset.useCondDB   = False
#process.ak5CaloL1Fastjet.useCondDB  = False
#
#process.ak5JPTL1Offset.useCondDB    = False
#process.ak5JPTL1Fastjet.useCondDB   = False
#
#process.ak5PFL1Offset.useCondDB     = False
#process.ak5PFL1Fastjet.useCondDB    = False
print "!!!===========================================================!!!"
print "!!!          dropping ak5L1JPTOffset                          !!!"
print "!!!          dropping ak5JPTL1Offset                          !!!"
print "!!!  this is not supposed to be necessary in a complete 44X   !!!"
print "!!!===========================================================!!!"
process.ak5JPTL2L3.correctors.remove('ak5L1JPTOffset')
process.ak5JPTL1L2L3.correctors.remove('ak5L1JPTOffset')
process.ak5JPTL1L2L3.correctors.remove('ak5JPTL1Offset')
process.ak5JPTL1FastL2L3.correctors.remove('ak5JPTL1Offset')
#process.ak5L1JPTOffset = cms.ESSource("LXXXCorrectionService",
#                                          useCondDB = cms.untracked.bool(False),
#                                          era = cms.string('Jec10V3'),
#                                          section = cms.string(''),
#                                          algorithm = cms.string('AK5JPT'),
#                                          level = cms.string('L1JPTOffset')
#                                      )

#undo what's pulled in by including Reconstruction_cff
#it relies on transient steps introduced in PF in 44X (back-fill)
process.pfPileUp.PFCandidates = cms.InputTag("particleFlow")
process.pfNoPileUp.bottomCollection = cms.InputTag("particleFlow") 

metJESCorAK5CaloJet.inputUncorJetsLabel = cms.string("ak5CaloJets")

# Input
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(
      'file:/nfs-3/userdata/cms2/cms2_validation/CMSSW_4_2_7_patch1/DoubleMuon/488D4763-00B1-E011-B110-BCAEC532971C.root'   # CMSSW_4_2_7
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
process.options                       = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )
process.source.noEventSort            = cms.untracked.bool( True )
process.MessageLogger.cerr.threshold  = ''

# Number of Events to Process
process.maxEvents                     = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Hypothesis cuts
process.hypDilepMaker.TightLepton_PtCut  = cms.double(20.0)
process.hypDilepMaker.LooseLepton_PtCut  = cms.double(10.0)
process.hypTrilepMaker.TightLepton_PtCut = cms.double(20.0)
process.hypTrilepMaker.LooseLepton_PtCut = cms.double(10.0)

# Event Maker
process.eventMaker.datasetName = cms.string("")
process.eventMaker.CMS2tag     = cms.string("")

