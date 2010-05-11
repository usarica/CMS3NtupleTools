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
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

process.load("RecoJets.Configuration.RecoJPTJets_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.GlobalTag.globaltag = "START3X_V26::All"

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000



#-------------------------------------------------
# PAT configuration
#-------------------------------------------------
#process.load("PhysicsTools.PatAlgos.patSequences_cff")
#process.patDefaultSequence = cms.Sequence(
#    process.patCandidates *
#    process.selectedPatCandidates
#)
#
##add muon isolation
#from PhysicsTools.PatAlgos.tools.muonTools import *
#addMuonUserIsolation.isolationTypes = ['All']
#addMuonUserIsolation.toolCode(process)
#
##change JetID tag
#from PhysicsTools.PatAlgos.tools.jetTools import *
#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
#run36xOn35xInput(process)
#addJetID( process, cms.InputTag('prunedUncorrectedCMS2Jets'), "antikt5" )
#switchJetCollection35X(process, 
#                    cms.InputTag('prunedUncorrectedCMS2Jets'),   
#                    doJTA            = True,            
#                    doBTagging       = True,            
#                    jetCorrLabel     = None,
#                    doType1MET       = True,
#                    genJetCollection = cms.InputTag("cms2antikt5GenJets"),
#                    doJetID          = True,
#                    jetIdLabel       = "cms2ak5"
#                    )

# add statement to prevent the PAT from using generator information
#from PhysicsTools.PatAlgos.tools.coreTools import *
#uncomment for data
#removeMCMatching(process, ['All'])

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
	'/store/mc/Spring10/QCD_Pt30/GEN-SIM-RECO/START3X_V26_S09-v1/0043/E09CE6E2-0947-DF11-90A9-001A4BA94900.root'
    #'file:/store/disk00/jribnik/Spring10_TTbarJets-madgraph_GEN-SIM-RECO_START3X_V26_S09-v1_0005_2AA58B20-AD46-DF11-9274-003048C69032.root'
    ),
)


#single lepton filter
process.EventSelectionJetPhotonFilt = cms.PSet(
    SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('pWithJetPhoton')
                )
        )

process.out = cms.OutputModule(
        "PoolOutputModule",
        process.EventSelectionJetPhotonFilt,
        verbose = cms.untracked.bool(True),
        dropMetaData = cms.untracked.string("NONE"),
        fileName = cms.untracked.string('ntuple.root')
)

process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))


# load event level configurations
process.load("CMS2.NtupleMaker.cms2CoreSequencesReduced_cff")
process.load('CMS2.NtupleMaker.jetphotonFilter_cfi')

# loosen thresholds on collections, change default jet cut
process.prunedUncorrectedCMS2Jets.usecorrectedCut = cms.bool(True)
process.scMaker.scEtMin = cms.double(5.0) #default, but maker sure
process.photonMaker.minEt = cms.double(5.0) #gev, min to keep--lower than default


#-------------------------------------------------
# process paths;
#-------------------------------------------------
process.cms2WithEverything             = cms.Sequence( process.cms2CoreSequence )


process.eventMaker.datasetName = cms.string("")
process.eventMaker.CMS2tag     = cms.string("")
#since filtering is done in the last step, there is no reason to remove these paths
#just comment out/remove an output which is not needed
process.pWithJetPhoton = cms.Path(process.cms2WithEverything * process.jetphotonFilter   )

#stuff to speed up I/O from castor
process.AdaptorConfig = cms.Service("AdaptorConfig",
                                  stats = cms.untracked.bool(True),
                                  enable = cms.untracked.bool(True),
                                  cacheHint = cms.untracked.string("lazy-download"),
                                  readHint = cms.untracked.string("auto-detect")
                                  )

process.source.noEventSort = cms.untracked.bool(True)


process.outpath         = cms.EndPath(process.out)


