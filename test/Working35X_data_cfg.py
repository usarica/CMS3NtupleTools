import FWCore.ParameterSet.Config as cms

process = cms.Process("CMS2")

from Configuration.EventContent.EventContent_cff import *

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.2 $'),
        annotation = cms.untracked.string('CMS2'),
        name = cms.untracked.string('CMS2 test configuration')
)

# load event level configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")
process.GlobalTag.globaltag = "GR09_R_35_V2B::All"

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

import sys
import os
import string


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 100


#-------------------------------------------------
# PAT configuration
#-------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.patDefaultSequence = cms.Sequence(
    process.patCandidates *
    process.selectedPatCandidates
)

#add muon isolation
from PhysicsTools.PatAlgos.tools.muonTools import *
addMuonUserIsolation.isolationTypes = ['All']
addMuonUserIsolation.toolCode(process)

#change JetID tag
from PhysicsTools.PatAlgos.tools.jetTools import *
#addJetID( process, cms.InputTag('prunedUncorrectedCMS2Jets'), "antikt5" )
switchJetCollection(process, 
                    cms.InputTag('prunedUncorrectedCMS2Jets'),   
                    doJTA            = True,            
                    doBTagging       = True,            
                    jetCorrLabel     = ('AK5','Calo'),  
                    doType1MET       = True,
                    genJetCollection = cms.InputTag("cms2antikt5GenJets"),
                    doJetID          = True,
#                    jetIdLabel       = "antikt5"
                    jetIdLabel       = "cms2ak5"
                    )

# add statement to prevent the PAT from using generator information
from PhysicsTools.PatAlgos.tools.coreTools import *
#uncomment for data
removeMCMatching(process, ['All'])

#-----------------------------------------------------------
# configure input data files and number of event to process
#-----------------------------------------------------------

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
'file:/data/tmp/kalavase/blah.root'
    ),
)




# load event level configurations
process.load("CMS2.NtupleMaker.cms2CoreSequences_cff")
process.load("CMS2.NtupleMaker.cms2GENSequence_cff")
process.load("CMS2.NtupleMaker.cms2PATSequence_cff")


# loosen thresholds on collections
process.hypDilepMaker.TightLepton_PtCut=cms.double(1)
process.hypDilepMaker.LooseLepton_PtCut=cms.double(1)

process.load("PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi")
#process.load("PhysicsTools.PatAlgos.recoLayer0.jetMETCorrections_cff")
#process.patJetCorrections = process.process.jetCorrFactors
#process.patJetCorrFactors = process.jetCorrFactors
#process.prunedUncorrectedCMS2Jets.inputUncorrectedJetCollection = cms.InputTag("antikt5CaloJets")
#process.pfJetMaker.pfJetsInputTag = cms.InputTag("antikt5PFJets")
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_7TeV_cff")
process.jetMaker.correctionTags   = cms.string("Summer09_7TeV_L2Relative_AK5Calo:Summer09_7TeV_L3Absolute_AK5Calo")
process.scjetMaker.correctionTags   = cms.string("Summer09_7TeV_L2Relative_SC5Calo:Summer09_7TeV_L3Absolute_SC5Calo")
#switchJECSet.newName = cms.string("Summer09_7TeV")
#switchJECSet.toolCode(process)
jetCorrFactors = getattr(process, 'patJetCorrFactors')
jetCorrFactors.corrSample = "Summer09_7TeV"

# don't forget 8e29
#process.hltMakerSequence += process.hlt8e29Maker
process.l1Maker.fillL1Particles = cms.untracked.bool(False)
##
process.load('CMS2.NtupleMaker.pixelDigiMaker_cfi')
process.load('CMS2.NtupleMaker.beamHaloSequence_cff')

process.RandomNumberGeneratorService.randomConeIsoMaker = cms.PSet( engineName = cms.untracked.string('HepJamesRandom'), 
        initialSeedSet = cms.untracked.vuint32(4126))

#-------------------------------------------------
# process paths;
#-------------------------------------------------
process.cms2WithEverythingExceptGEN  = cms.Sequence(   process.coreCMS2Sequence
                                                     * process.cms2beamHaloSequence
                                                     * process.pixelDigiMaker
                                                     * process.patDefaultSequence * process.cms2PATSequence)

#since filtering is done in the last step, there is no reason to remove these paths
#just comment out/remove an output which is not needed
#process.pWithRecoLepton      = cms.Path(process.cms2WithEverythingExceptGEN * process.aSkimFilter )
process.pNoFilter = cms.Path(process.cms2WithEverythingExceptGEN)
process.eventMaker.datasetName = cms.string("/MinimumBias/BeamCommissioning09-RecoTracks-Mar3rdSkim_v2/RAW-RECO")
process.eventMaker.CMS2tag     = cms.string("blah")

process.EventSelectionSingleFilt = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pWithRecoLepton')
    )
)
process.outRecoOrGenSingleFilt_CMS2 = cms.OutputModule(
        "PoolOutputModule",
        process.EventSelectionSingleFilt,
        verbose = cms.untracked.bool(True),
        dropMetaData = cms.untracked.string("NONE"),
        fileName = cms.untracked.string('ntuple.root')
)
process.outRecoOrGenSingleFilt_CMS2.outputCommands = cms.untracked.vstring( 'drop *' )
process.outRecoOrGenSingleFilt_CMS2.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))

process.outNoFilter_CMS2 = cms.OutputModule(
        "PoolOutputModule",
        verbose = cms.untracked.bool(True),
        dropMetaData = cms.untracked.string("NONE"),
        fileName = cms.untracked.string('ntuple.root')
)
process.outNoFilter_CMS2.outputCommands = cms.untracked.vstring( 'drop *' )
process.outNoFilter_CMS2.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))

#process.outpath         = cms.EndPath(process.outRecoOrGenSingleFilt_CMS2)
process.outpath         = cms.EndPath(process.outNoFilter_CMS2)
