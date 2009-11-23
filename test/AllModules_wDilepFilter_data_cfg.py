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
process.GlobalTag.globaltag = "MC_31X_V3::All"

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)



process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1



#-------------------------------------------------
# PAT configuration
#-------------------------------------------------

process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.patDefaultSequence = cms.Sequence(
    process.allLayer1Objects *
    process.selectedLayer1Objects
)


#change JetID tag
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process, 
                    cms.InputTag('prunedUncorrectedCMS2Jets'),   
                    doJTA            = True,            
                    doBTagging       = True,            
                    jetCorrLabel     = ('AK5','Calo'),  
                    doType1MET       = True,
                    genJetCollection = cms.InputTag("cms2antikt5GenJets"),
                    doJetID          = False,
                    jetIdLabel       = "ak5"
                    )

# add statement to prevent the PAT from using generator information
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, 'All')

#-----------------------------------------------------------
# configure input data files and number of event to process
#-----------------------------------------------------------

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
        'file:TTbar_333_RelVal_1.root',
        'file:TTbar_333_RelVal_2.root'
        #'file:TTbar_333_RelVal_3.root'
    ),
    #won't need this when running on actual data 
    inputCommands = cms.untracked.vstring(
    'keep *',
    'drop *_antikt5GenJets_*_*',
    'drop *_g4SimHits_*_*',
    'drop *_generator_*_*',
    'drop *_genMetCalo_*_*',
    'drop *_genMetCaloAndNonPrompt_*_*',
    'drop *_genMetTrue_*_*',
    'drop *_genParticles_*_*',
    'drop *_iterativeCone5GenJets_*_*',
    'drop *_kt4GenJets_*_*',
    'drop *_kt6GenJets_*_*',
    'drop *_simMuonCSCDigis_*_*',
    'drop *_simMuonDTDigis_*_*',
    'drop *_simMuonRPCDigis_*_*',
    'drop *_sisCone5GenJets_*_*',
    'drop *_sisCone7GenJets_*_*'
    )
)

## define event selection
process.EventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)

process.out_CMS2 = cms.OutputModule(
        "PoolOutputModule",
    process.EventSelection,
    verbose = cms.untracked.bool(True),
    dropMetaData = cms.untracked.string("NONE"),
    fileName = cms.untracked.string('ntuple_AllModules_wDilepFilter_data.root')
)

process.out_CMS2.outputCommands = cms.untracked.vstring( 'drop *' )
process.out_CMS2.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))


# load event level configurations
process.load("CMS2.NtupleMaker.cms2CoreSequences_cff")
process.load("CMS2.NtupleMaker.cms2PATSequence_cff")
process.load("CMS2.NtupleMaker.hypFilter_cfi")

#-------------------------------------------------
# process paths;
#-------------------------------------------------
process.all             = cms.Sequence( process.coreCMS2Sequence * process.patDefaultSequence * process.cms2PATSequence)

process.p              = cms.Path(process.all * process.hypFilter)

process.outpath       = cms.EndPath(process.out_CMS2)

