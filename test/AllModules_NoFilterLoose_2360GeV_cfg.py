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
process.GlobalTag.globaltag = "MC_31X_V3::All"

process.options = cms.untracked.PSet(
   Rethrow = cms.untracked.vstring('ProductNotFound')
)



process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1
#use 900 GeV JEC
process.load("JetMETCorrections.Configuration.L2L3Corrections_2360GeV_cff")


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
#        '/store/data/BeamCommissioning09/MinimumBias/RECO/rereco_FIRSTCOLL_v1/0083/FE5EDBBC-7DD9-DE11-9589-001A92971B64.root'
'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec9thReReco-v1/0002/1E41E0CE-C2E4-DE11-B1CE-0026189438D7.root'
    ),
    eventsToProcess = cms.untracked.VEventRange("123596:15707761-123596:MAX"),
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
    fileName = cms.untracked.string('ntuple_2360.root')
)

process.out_CMS2.outputCommands = cms.untracked.vstring( 'drop *' )
process.out_CMS2.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))


# load event level configurations
process.load("CMS2.NtupleMaker.cms2CoreSequences_cff")



# loosen thresholds on collections
process.scMaker.scEtMin = cms.double(0.0)
process.trkJetMaker.trkJetPtCut = cms.double(0.0)
process.prunedUncorrectedCMS2Jets.uncorrectedJetPtCut = cms.double(0.0)
process.jetMaker.correctionLevels = cms.string("L2:L3")
process.jetMaker.correctionTags   = cms.string("2360GeV_L2Relative_AK5Calo:2360GeV_L3Absolute_AK5Calo")
process.hypDilepMaker.TightLepton_PtCut = cms.double(0.0)
process.hypDilepMaker.LooseLepton_PtCut = cms.double(0.0)
process.hypTrilepMaker.TightLepton_PtCut = cms.double(0.0)
process.hypTrilepMaker.LooseLepton_PtCut = cms.double(0.0)
process.hypQuadlepMaker.TightLepton_PtCut = cms.double(0.0)
process.hypQuadlepMaker.LooseLepton_PtCut = cms.double(0.0)

#
# set up random seeds
#

process.RandomNumberGeneratorService.randomConeIsoMaker = cms.PSet( engineName = cms.untracked.string('HepJamesRandom'), 
	initialSeedSet = cms.untracked.vuint32(4126))

#-------------------------------------------------
# process paths;
#-------------------------------------------------
process.p             = cms.Path( process.coreCMS2Sequence)

process.outpath       = cms.EndPath(process.out_CMS2)

