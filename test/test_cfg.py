# The following comments couldn't be translated into the new config version:

# -*-sh-*-

import FWCore.ParameterSet.Config as cms

process = cms.Process("CMS2")

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.20 $'),
        annotation = cms.untracked.string('CMS2'),
        name = cms.untracked.string('CMS2 test configuration')
)

## configure geometry
process.load("Configuration.StandardSequences.Geometry_cff")
## configure B field
process.load("Configuration.StandardSequences.MagneticField_cff")
## configure conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "IDEAL_V9::All"

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''

#----------------------------------------------------
#CMS2 includes
#----------------------------------------------------
process.load("CMS2.NtupleMaker.beamSpotMaker_cfi")
process.load("CMS2.NtupleMaker.eventMaker_cfi")
process.load("CMS2.NtupleMaker.jetMaker_cfi")
process.load("CMS2.NtupleMaker.jetToMuAssMaker_cfi")
process.load("CMS2.NtupleMaker.jetToElAssMaker_cfi")
process.load("CMS2.NtupleMaker.patJetMaker_cfi")
process.load("CMS2.NtupleMaker.electronMaker_cfi")
process.load("CMS2.NtupleMaker.patElectronMaker_cfi")
process.load("CMS2.NtupleMaker.hypTrilepMaker_cfi")
#process.load("CMS2.NtupleMaker.metMaker_cfi")
process.load("CMS2.NtupleMaker.metSequence_cff")
process.load("CMS2.NtupleMaker.genMaker_cfi")
process.load("CMS2.NtupleMaker.candToGenAssMaker_cfi")
process.load("CMS2.NtupleMaker.muonMaker_cfi")
process.load("CMS2.NtupleMaker.trackMaker_cfi")
process.load("CMS2.NtupleMaker.patMETMaker_cfi")
process.load("CMS2.NtupleMaker.patMuonMaker_cfi")
process.load("CMS2.NtupleMaker.muToElsAssMaker_cfi")
process.load("CMS2.NtupleMaker.muToJetAssMaker_cfi")
process.load("CMS2.NtupleMaker.muToTrackAssMaker_cfi")
process.load("CMS2.NtupleMaker.elToTrackAssMaker_cfi")
process.load("CMS2.NtupleMaker.elToMuAssMaker_cfi")
process.load("CMS2.NtupleMaker.trackToMuonAssMaker_cfi")
process.load("CMS2.NtupleMaker.trackToElsAssMaker_cfi")
process.load("CMS2.NtupleMaker.hypDilepMaker_cfi")
process.load("CMS2.NtupleMaker.hypQuadlepMaker_cfi")
process.load("CMS2.NtupleMaker.triggerEventMaker_cfi")
process.load("CMS2.NtupleMaker.l1DigiMaker_cfi")
process.load("CMS2.NtupleMaker.theFilter_cfi")
process.load("CMS2.NtupleMaker.elCaloIsoSequence_cff")
process.load("CMS2.NtupleMaker.genJetMaker_cfi")
process.load("CMS2.NtupleMaker.conversionMaker_cfi")
process.load("CMS2.NtupleMaker.genJetMaker_cfi")
process.load("CMS2.NtupleMaker.trkMuonFilter_cfi")
process.load("CMS2.NtupleMaker.trkJetMaker_cfi")
process.load("CMS2.NtupleMaker.tcmetMaker_cfi")
#process.Timing = cms.Service("Timing")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

##source 
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring('file:///uscms_data/d1/slava77/tauola-16AAC418-218A-DD11-AC33-001F2908F0E4.root')
)


#-------------------------------------------------
# patTuple configuration
#-------------------------------------------------

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopObjectProducers.patTuple_cff")

## necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)


#-------------------------------------------------
# JEC
#-------------------------------------------------
#############   Define the L2 correction service #####
process.L2RelativeJetCorrector = cms.ESSource("L2RelativeCorrectionService", 
    tagName = cms.string('Summer08_L2Relative_SC5Calo'),
    label = cms.string('L2RelativeJetCorrector')
)
#############   Define the L3 correction service #####
process.L3AbsoluteJetCorrector = cms.ESSource("L3AbsoluteCorrectionService", 
    tagName = cms.string('Summer08_L3Absolute_SC5Calo'),
    label = cms.string('L3AbsoluteJetCorrector')
)
#############   Define the EMF correction service ####
process.L4JetCorrector = cms.ESSource("L4EMFCorrectionService", 
    tagName = cms.string('CMSSW_152_L4EMF'),
    label = cms.string('L4EMFJetCorrector')
)
#############   Define the chain corrector service - L2L3 ###
process.L2L3JetCorrector = cms.ESSource("JetCorrectionServiceChain",  
    correctors = cms.vstring('L2RelativeJetCorrector','L3AbsoluteJetCorrector'),
    label = cms.string('L2L3JetCorrector')
)
#############   Define the chain corrector service - L2L3L4 ###
process.L2L3L4JetCorrector = cms.ESSource("JetCorrectionServiceChain",  
    correctors = cms.vstring('L4EMFJetCorrector'),
    label = cms.string('L2L3L4JetCorrector')
)
#############   Define the chain corrector module - L2L3L4####
process.L2L3CorJet = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("sisCone5CaloJets"),
    correctors = cms.vstring('L2L3JetCorrector')
)
#############   Define the chain corrector module - L2L3L4####
process.L2L3L4CorJet = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("L2L3CorJet"), ##don't want to rerun L2L3 for L2L3L4 -> take L2L3 as input
    correctors = cms.vstring('L2L3L4JetCorrector')
)
# set the record's IOV. Must be defined once. Choose ANY correction service. #
process.prefer("L2L3L4JetCorrector") 


## change jet collection
from PhysicsTools.PatAlgos.tools.jetTools import *

switchJetCollection(process, 
                    'sisCone5CaloJets',             # jet collection; must be already in the event when patLayer0 sequence is executed
                    layers       = [0,1],           # if you're not running patLayer1, set 'layers=[0]' 
                    runCleaner   = "CaloJet",       # =None if not to clean
                    doJTA        = True,            # run jet-track association & JetCharge
                    doBTagging   = True,            # run b-tagging
                    jetCorrLabel = ('SC5', 'Calo'), # example jet correction name; set to None for no JEC
                    doType1MET   = True             # recompute Type1 MET using these jets
                    )

#-------------------------------------------------
# load tcMET producer
#-------------------------------------------------

process.load("RecoMET.METProducers.TCMET_cfi")

#-------------------------------------------------
# process output; first the event selection is
# defined: only those events that have passed the
# full production path are selected and written
# to file; the event content has been defined
# above
#-------------------------------------------------

## define event selection
process.EventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)


##default Jurassic Isolation values
##should be the same as in the PAT, but just being safe
process.allLayer1Electrons.isolation.ecal.vetos = cms.vstring(
    'EcalBarrel:0.045',
    'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
    'EcalBarrel:ThresholdFromTransverse(0.08)',
    'EcalEndcaps:ThresholdFromTransverse(0.3)',
    'EcalEndcaps:0.070',
    'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)'
)

process.allLayer0Electrons.isolation.ecal.vetos = cms.vstring(
    'EcalBarrel:0.045',
    'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
    'EcalBarrel:ThresholdFromTransverse(0.08)',
    'EcalEndcaps:ThresholdFromTransverse(0.3)',
    'EcalEndcaps:0.070',
    'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)'
)


process.l1DigiMaker = cms.EDFilter("L1DigiMaker")

process.CMS2 = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'drop *',
    ##'keep recoGenParticles_genParticles_*_*',
    ##'keep *_genEventScale_*_*',
    ##'keep *_genEventWeight_*_*',
    ##'keep *_genEventPdfInfo_*_*',
    ##'keep edmTriggerResults_TriggerResults_*_HLT', 
    ##'keep *_hltTriggerSummaryAOD_*_*',
    ##'keep *_offlineBeamSpot_*_*',
    ##'keep *_offlinePrimaryVertices_*_*',
    ##'keep recoTracks_generalTracks_*_*', 
    ##'keep *_towerMaker_*_*',
    ##'keep *_selectedLayer1Photons_*_*', 
    'keep *_selectedLayer1Electrons_*_*', 
    ##'keep *_selectedLayer1Muons_*_*', 
    ##'keep *_selectedLayer1Taus_*_*', 
    ##'keep *_selectedLayer1Jets_*_*', 
    ##'keep *_selectedLayer1METs_*_*',
    ##'keep patPFParticles_*_*_*',
    ##'keep *_selectedLayer1Hemispheres_*_*',
    'keep *_*Maker_*_CMS2',
    'keep recoCaloMET*_*_*_*'
  )
)

## configure output module
process.out = cms.OutputModule("PoolOutputModule",
    process.EventSelection,
    process.CMS2,
    ##process.patTupleEventContent,
    verbose = cms.untracked.bool(True),
    dropMetaDataForDroppedData = cms.untracked.bool(True),                           
    fileName = cms.untracked.string('test_1.root')
)


#-------------------------------------------------
# process paths;
#-------------------------------------------------

#process.triggerEventMaker = cms.EDProducer("TriggerEventMaker")

process.MetCorrection = cms.Sequence(process.tcMet*process.tcmetMaker)
process.JetCorrection = cms.Sequence(process.L2L3CorJet*process.L2L3L4CorJet)
process.makers = cms.Sequence(process.beamSpotMaker*process.muonMaker*process.electronMaker*process.jetMaker*process.trackMaker)
process.patmakers = cms.Sequence(process.patMuonMaker*process.patElectronMaker*process.patJetMaker*process.patMETMaker)
process.assmakers = cms.Sequence(process.jetToMuAssMaker*process.jetToElAssMaker*process.muToElsAssMaker*process.candToGenAssMaker*process.muToJetAssMaker*process.muToTrackAssMaker*process.elToTrackAssMaker*process.elToMuAssMaker*process.trackToMuonAssMaker*process.trackToElsAssMaker)
process.trigprimmakers = cms.Sequence(process.l1DigiMaker*process.triggerEventMaker)
#process.generalmakers = cms.Sequence(process.eventMaker*process.metMaker*process.genMaker*process.genjetmaker)
process.generalmakers = cms.Sequence(process.eventMaker*process.metCorSequence*process.genMaker*process.genjetmaker)
process.hypmaker = cms.Sequence(process.hypTrilepMaker*process.hypDilepMaker*process.hypQuadlepMaker)
process.othermakers = cms.Sequence(process.elCaloIsoSequence*process.conversionMaker)
process.cms2 = cms.Sequence(process.generalmakers*process.trigprimmakers*process.makers*process.patmakers*process.assmakers*process.hypmaker*process.genjetmaker*process.trkmuonfilter*process.trkjetmaker*process.othermakers)

##includes filter
#process.p = cms.Path(process.MetCorrection*process.JetCorrection*process.patTuple*process.cms2*process.theFilter)

##no filter
process.p = cms.Path(process.MetCorrection*process.JetCorrection*process.patTuple*process.cms2)

##output
process.outpath = cms.EndPath(process.out)


