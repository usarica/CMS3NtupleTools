import FWCore.ParameterSet.Config as cms

process = cms.Process("CMS2")

from Configuration.EventContent.EventContent_cff import *

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.1 $'),
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

process.load("RecoJets.Configuration.RecoJPTJets_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('Configuration/EventContent/EventContent_cff')

#process.GlobalTag.globaltag = "GR09_R_35_V2B::All"
# changed at 16:45 Apr 7
process.GlobalTag.globaltag = "GR10_E_V5::All"

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
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
run36xOn35xInput(process)
addJetID( process, cms.InputTag('prunedUncorrectedCMS2Jets'), "antikt5" )
switchJetCollection(process, 
                    cms.InputTag('prunedUncorrectedCMS2Jets'),   
                    doJTA            = True,            
                    doBTagging       = True,            
                    jetCorrLabel     = ('AK5','Calo'),  
                    doType1MET       = True,
                    genJetCollection = cms.InputTag("cms2antikt5GenJets"),
                    doJetID          = True,
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
'rfio:///castor/cern.ch/cms/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/0065C919-F53B-DF11-8BF5-001D09F29146.root'
#'file:/home/users/jribnik/devel/minbias352.root'
    ),
)




# load event level configurations
process.load("CMS2.NtupleMaker.cms2CoreSequences_cff")
process.load("CMS2.NtupleMaker.cms2PATSequence_cff")
process.load("CMS2.NtupleMaker.caloTowerSequence_cff")

process.hltMaker.processName = cms.untracked.string("HLT")
process.hltMakerSequence = cms.Sequence(process.hltMaker)

# loosen thresholds on collections
process.hypDilepMaker.TightLepton_PtCut = cms.double(0)
process.hypDilepMaker.LooseLepton_PtCut = cms.double(0)
process.hypDilepMaker.useSTAMuon = cms.bool(True)
process.eventMaker.isData      = cms.bool(True)
process.eventMaker.datasetName = cms.string("/MinimumBias/BeamCommissioning09-RecoTracks-Mar3rdSkim_v2/RAW-RECO")
process.eventMaker.CMS2tag     = cms.string("blah")
process.l1Maker.fillL1Particles = cms.untracked.bool(True)
process.l1Maker.l1ParticlesProcessName = cms.untracked.string('EXPRESS')

process.RandomNumberGeneratorService.randomConeIsoMaker = cms.PSet( engineName = cms.untracked.string('HepJamesRandom'), 
        initialSeedSet = cms.untracked.vuint32(4126))

#-------------------------------------------------
# process paths;
#-------------------------------------------------
process.cms2WithEverythingExceptGEN  = cms.Sequence(   process.cms2CoreSequence
                                                     * process.cms2CaloTowerSequence  
                                                     * process.patDefaultSequence
                                                     * process.cms2PATSequence
                                                     * process.randomConeIsoMaker  )

#since filtering is done in the last step, there is no reason to remove these paths
#just comment out/remove an output which is not needed
#process.pWithRecoLepton      = cms.Path(process.cms2WithEverythingExceptGEN * process.aSkimFilter )
process.pNoFilter = cms.Path(process.cms2WithEverythingExceptGEN)

process.EventSelectionSingleFilt = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pWithRecoLepton')
    )
)

process.outNoFilter_CMS2 = cms.OutputModule(
        "PoolOutputModule",
        verbose = cms.untracked.bool(True),
        dropMetaData = cms.untracked.string("NONE"),
        fileName = cms.untracked.string('ntuple.root')
)
process.outNoFilter_CMS2.outputCommands = cms.untracked.vstring( 'drop *' )
process.outNoFilter_CMS2.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))

process.outpath         = cms.EndPath(process.outNoFilter_CMS2)


process.AdaptorConfig = cms.Service("AdaptorConfig",
                                    stats = cms.untracked.bool(True),
                                    enable = cms.untracked.bool(True),
                                    cacheHint = cms.untracked.string("lazy-download"),
                                    readHint = cms.untracked.string("auto-detect")
                                    )

process.source.fileNames = cms.untracked.vstring(os.environ['INPUT_FILE'])
process.source.skipEvents = cms.untracked.uint32(string.atoi(os.environ['SKIP_EVENTS']))
process.source.noEventSort = cms.untracked.bool(True)
process.outNoFilter_CMS2.fileName = cms.untracked.string(os.environ['OUTPUT_FILE'])
