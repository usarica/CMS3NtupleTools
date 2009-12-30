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
	'/store/mc/Summer09/ppMuXLoose/GEN-SIM-RECO/STARTUP3X_V8F_2360GeV-v1/0088/06C46427-75E0-DE11-879B-002618943973.root'	
	#'/store/data/BeamCommissioning09/MinimumBias/RECO/rereco_FIRSTCOLL_v1/0083/FE5EDBBC-7DD9-DE11-9589-001A92971B64.root'
	#'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec9thReReco-v1/0002/1E41E0CE-C2E4-DE11-B1CE-0026189438D7.root'
    ),
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
process.load("CMS2.NtupleMaker.cms2CoreSequences_cff")
process.load("CMS2.NtupleMaker.cms2GENSequence_cff")


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
process.p             = cms.Path( process.coreCMS2Sequence * process.cms2GENSequence )

process.outpath       = cms.EndPath(process.out_CMS2)

