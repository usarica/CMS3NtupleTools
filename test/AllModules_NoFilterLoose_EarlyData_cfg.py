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


#-----------------------------------------------------------
# configure input data files and number of event to process
#-----------------------------------------------------------

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
    'file:BeamCommissioning09_MinimumBias_0083_022E364E-28D9-DE11-985F-0026189438B4.root'
    #'file:/home/users/kalavase/work/MinBiasData_900GeVRun/V03-00-07/CMSSW_3_3_3/src/CMS2/NtupleMaker/test/TTbar_333_RelVal_1.root',
     #   'file:/home/users/kalavase/work/MinBiasData_900GeVRun/V03-00-07/CMSSW_3_3_3/src/CMS2/NtupleMaker/test/TTbar_333_RelVal_2.root'
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
    'drop *_sisCone7GenJets_*_*',
    'drop *_genMetIC5GenJets_*_*',
    'drop *_genMetIC7GenJets_*_*'
    )
)


process.out_CMS2 = cms.OutputModule(
        "PoolOutputModule",
    verbose = cms.untracked.bool(True),
    dropMetaData = cms.untracked.string("NONE"),
    fileName = cms.untracked.string('ntuple.root')
)

process.out_CMS2.outputCommands = cms.untracked.vstring( 'drop *' )
process.out_CMS2.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))


# load event level configurations
process.load("CMS2.NtupleMaker.cms2CoreSequences_EarlyData_cff")

# loosen thresholds on collections
process.scMaker.scEtMin = cms.double(0.0)
process.trkJetMaker.trkJetPtCut = cms.double(0.0)
process.prunedUncorrectedCMS2Jets.uncorrectedJetPtCut = cms.double(0.0)

#-------------------------------------------------
# process paths;
#-------------------------------------------------
process.p             = cms.Path( process.coreCMS2Sequence_EarlyData)

process.outpath       = cms.EndPath(process.out_CMS2)

