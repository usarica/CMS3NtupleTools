import FWCore.ParameterSet.Config as cms

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#VarParsing_Documentation
# Allow command line options like
#     cmsRun main_pset.py data=True prompt=True   # prompt data
#     cmsRun main_pset.py data=False               # MC
#     cmsRun main_pset.py fastsim=True             # fastsim
import FWCore.ParameterSet.VarParsing as VarParsing
opts = VarParsing.VarParsing('python')
vpbool = VarParsing.VarParsing.varType.bool
opts.register('data'    , False  , mytype=vpbool)
opts.register('prompt'  , False  , mytype=vpbool)
opts.register('fastsim' , False , mytype=vpbool)
opts.register('relval'  , False , mytype=vpbool)
opts.register('triginfo'  , False , mytype=vpbool)
opts.parseArguments()
# be smart. if fastsim, it's obviously MC
# if it's MC, it's obviously not prompt
if opts.fastsim: opts.data = False
if not opts.data: opts.prompt = False
print """PSet is assuming:
   data? {}
   prompt? {}
   fastsim? {}
   relval? {}
   triginfo? {}
""".format(bool(opts.data), bool(opts.prompt), bool(opts.fastsim), bool(opts.relval), bool(opts.triginfo))

import CMS3.NtupleMaker.configProcessName as configProcessName
configProcessName.name="PAT"
if opts.data and opts.prompt:
    configProcessName.name="RECO"

configProcessName.name2="RECO"

if opts.relval:
    if opts.data:
        configProcessName.name="reRECO"
        configProcessName.name2="reRECO"
    else:
        configProcessName.name="RECO"
        configProcessName.name2="RECO"

if opts.fastsim:
    configProcessName.fastSimName="HLT"
    configProcessName.name2=configProcessName.fastSimName
configProcessName.isFastSim=opts.fastsim

# CMS3
process = cms.Process("CMS3")

# Version Control For Python Configuration Files
process.configurationMetadata = cms.untracked.PSet(
        version    = cms.untracked.string('$Revision: 1.11 $'),
        annotation = cms.untracked.string('CMS3'),
        name       = cms.untracked.string('CMS3 test configuration')
)

from Configuration.EventContent.EventContent_cff   import *

# load event level configurations
process.load('Configuration/EventContent/EventContent_cff')
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

# services
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.GlobalTag.globaltag = "80X_mcRun2_asymptotic_2016_miniAODv2_v0" #80X
#process.GlobalTag.globaltag = "91X_upgrade2017_realistic_v5" #MC
#process.GlobalTag.globaltag = "91X_dataRun2_relval_v6" #data
process.GlobalTag.globaltag = "94X_mc2017_realistic_v14"
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold  = ''
process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')
# process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True),SkipEvent = cms.untracked.vstring('ProductNotFound') )
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(False) )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName     = cms.untracked.string('ntuple.root'),
                               # fileName     = cms.untracked.string('ntuple2.root'),
                               dropMetaData = cms.untracked.string("ALL"),
                               # basketSize = cms.untracked.int32(16384*150)
                               # basketSize = cms.untracked.int32(16384*10)
                               basketSize = cms.untracked.int32(16384*23)
)

########    override the GT for MC     ########
### ESPrefer for L1TGlobalPrescalesVetosRcd ###
if not opts.data:
    process.load("CondCore.CondDB.CondDB_cfi")
    process.CondDB.connect = "frontier://FrontierProd/CMS_CONDITIONS"
    process.l1tPS = cms.ESSource("PoolDBESSource",
        process.CondDB,
        toGet = cms.VPSet(
            cms.PSet(
            record = cms.string("L1TGlobalPrescalesVetosRcd"),
            tag = cms.string("L1TGlobalPrescalesVetos_passThrough_mc")
            )
        )
    )
    process.es_prefer_l1tPS = cms.ESPrefer("PoolDBESSource", "l1tPS")
# tag from https://cms-conddb.cern.ch/cmsDbBrowser/list/Prod/gts/92X_upgrade2017_realistic_v2

#load cff and third party tools
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from JetMETCorrections.Configuration.CorrectedJetProducersDefault_cff import *
from JetMETCorrections.Configuration.CorrectedJetProducers_cff import *
from JetMETCorrections.Configuration.CorrectedJetProducersAllAlgos_cff import *
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
# from RecoJets.JetProducers.fixedGridRhoProducerFastjet_cfi import *
# process.fixedGridRhoFastjetAll = fixedGridRhoFastjetAll.clone(pfCandidatesTag = 'packedPFCandidates')

#Electron Identification for PHYS 14
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *  
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons',"",configProcessName.name)
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons',"",configProcessName.name)
process.egmGsfElectronIDSequence = cms.Sequence(process.electronMVAVariableHelper * process.electronMVAValueMapProducer * process.egmGsfElectronIDs)
my_id_modules = [
        # 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
        # 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
        # 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
        # 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
                 ]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# Load Ntuple producer cff
process.load("CMS3.NtupleMaker.cms3CoreSequences_cff")
if not opts.data: process.load("CMS3.NtupleMaker.cms3GENSequence_cff")
process.load("CMS3.NtupleMaker.cms3PFSequence_cff")
process.eventMaker.isData                        = cms.bool(opts.data)
    
# if do_deepbtag:
#     from PhysicsTools.PatAlgos.tools.jetTools import *
#     deep_discriminators = ["pfDeepCSVJetTags:probudsg", "pfDeepCSVJetTags:probb", "pfDeepCSVJetTags:probc", "pfDeepCSVJetTags:probbb", "pfDeepCSVJetTags:probcc" ]
#     updateJetCollection(
#         process,
#         jetSource = cms.InputTag('slimmedJets'),
#        jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
#         btagDiscriminators = deep_discriminators
#     )
#     updateJetCollection(
#         process,
#         labelName = 'Puppi',
#         jetSource = cms.InputTag('slimmedJetsPuppi'),
#        jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
#         btagDiscriminators = deep_discriminators
#     )

    # Needed for the above updateJetCollection() calls
    # process.pfJetMaker.pfJetsInputTag = cms.InputTag('selectedUpdatedPatJets')
    # process.pfJetPUPPIMaker.pfJetsInputTag = cms.InputTag('selectedUpdatedPatJetsPuppi')

# Hypothesis cuts
process.hypDilepMaker.TightLepton_PtCut  = cms.double(10.0)
process.hypDilepMaker.LooseLepton_PtCut  = cms.double(10.0)

#Options for Input
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'file:/home/users/dpgilber/2017/12CEC8EA-0743-E811-BE6A-0CC47A7C3430.root' # From: /DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
                                )
)
process.source.noEventSort = cms.untracked.bool( True )

#Max Events
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


#Run corrected MET maker

#configurable options =======================================================================
usePrivateSQlite=False #use external JECs (sqlite file)
applyResiduals=opts.data #application of residual corrections. Have to be set to True once the 13 TeV residual corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or for analysis checks and developments (not the official recommendation!).
#===================================================================

if usePrivateSQlite:
    from CondCore.DBCommon.CondDBSetup_cfi import *
    import os
    era="Summer15_25nsV5_MC"
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string( "sqlite_file:"+era+".db" ),
                               toGet =  cms.VPSet(
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                label= cms.untracked.string("AK4PF")
                ),
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                label= cms.untracked.string("AK4PFchs")
                ),
            )
                               )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

### =================================================================================
#jets are rebuilt from those candidates by the tools, no need to do anything else
### =================================================================================

process.outpath = cms.EndPath(process.out)
process.out.outputCommands = cms.untracked.vstring( 'drop *' )

if not opts.data:
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    #default configuration for miniAOD reprocessing, change the isData flag to run on data
    #for a full met computation, remove the pfCandColl input
    runMetCorAndUncFromMiniAOD(process,
                               isData=opts.data,
                               )

process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS3*'))

### -------------------------------------------------------------------
### the lines below remove the L2L3 residual corrections when processing data
### -------------------------------------------------------------------
if not applyResiduals:
    process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
### ------------------------------------------------------------------

# end Run corrected MET maker

# Extra trigger information (matching)
if opts.triginfo:
    process.load("CMS3.NtupleMaker.muToTrigAssMaker_cfi")
    process.load("CMS3.NtupleMaker.elToTrigAssMaker_cfi")
    if opts.data and opts.prompt:
        # process.muToTrigAssMaker.processName = cms.untracked.string("RECO")
        # process.elToTrigAssMaker.processName = cms.untracked.string("RECO")
        process.muToTrigAssMaker.triggerObjectsName = cms.untracked.string("slimmedPatTrigger")
        process.elToTrigAssMaker.triggerObjectsName = cms.untracked.string("slimmedPatTrigger")
        process.hltMaker.triggerObjectsName = cms.untracked.string("slimmedPatTrigger")
    process.hltMaker.fillTriggerObjects = cms.untracked.bool(True)

# python -c "from PhysicsTools.NanoAOD.electrons_cff import isoForEle; print 'process.isoForEle = {}'.format(repr(isoForEle))"
process.isoForEle = cms.EDProducer("EleIsoValueMapProducer",
    EAFile_MiniIso = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt'),
    EAFile_PFIso = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
    relative = cms.bool(False),
    rho_MiniIso = cms.InputTag("fixedGridRhoFastjetAll"),
    rho_PFIso = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedElectrons")
)
# python -c "from PhysicsTools.NanoAOD.muons_cff import isoForMu; print 'process.isoForMu = {}'.format(repr(isoForMu))"
process.isoForMu = cms.EDProducer("MuonIsoValueMapProducer",
    EAFile_MiniIso = cms.FileInPath('PhysicsTools/NanoAOD/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
    relative = cms.bool(False),
    rho_MiniIso = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedMuons")
)

process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName=cms.string('TransientTrackBuilder')
)

if opts.data:
    process.p = cms.Path( 
        process.metFilterMaker *
        process.egmGsfElectronIDSequence *     
        process.vertexMaker *
        process.secondaryVertexMaker *
        process.eventMaker *
        process.pfCandidateMaker *
        process.isoTrackMaker *
        process.isoForEle * 
        process.isoForMu *
        process.electronMaker *
        process.muonMaker *
        process.pfJetMaker *
        process.pfJetPUPPIMaker *
        process.subJetMaker *
        process.pfmetMaker *
        process.pfmetpuppiMaker *
        process.hltMakerSequence *
        process.miniAODrhoSequence *
        process.pftauMaker *
        process.photonMaker *
        # process.muToTrigAssMaker * # Note these are hacked in below
        # process.elToTrigAssMaker * # Note these are hacked in below
        # process.genMaker *
        # process.genJetMaker *
        # process.candToGenAssMaker * # requires electronMaker, muonMaker, pfJetMaker, photonMaker
        # process.pdfinfoMaker *
        # process.puSummaryInfoMaker *
        process.hypDilepMaker
    )
else:
    process.p = cms.Path( 
        process.metFilterMaker *
        process.egmGsfElectronIDSequence *     
        process.vertexMaker *
        process.secondaryVertexMaker *
        process.eventMaker *
        process.pfCandidateMaker *
        process.isoTrackMaker *
        process.isoForEle * 
        process.isoForMu *
        process.electronMaker *
        process.muonMaker *
        process.pfJetMaker *
        process.pfJetPUPPIMaker *
        process.subJetMaker *
        process.pfmetMaker *
        process.pfmetpuppiMaker *
        process.hltMakerSequence *
        process.miniAODrhoSequence *
        process.pftauMaker *
        process.photonMaker *
        process.genMaker *
        process.genJetMaker *
        process.candToGenAssMaker * # requires electronMaker, muonMaker, pfJetMaker, photonMaker
        process.pdfinfoMaker *
        process.puSummaryInfoMaker *
        process.hypDilepMaker
    )
if opts.triginfo:
    # Now insert the xToTrigAssMakers into the path
    # Hooray for hacky operations
    process.p.insert(process.p.index(process.photonMaker)+1,process.muToTrigAssMaker)
    process.p.insert(process.p.index(process.photonMaker)+1,process.elToTrigAssMaker)


process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.Timing = cms.Service("Timing",
        summaryOnly = cms.untracked.bool(True)
        )


# for use with Valgrind. After enabling, can do
# $ valgrind --leak-check=yes  cmsRun main_pset.py >& log.txt
# $ valgrindMemcheckParser.pl --preset=prod,-prod1+ log.txt  > blah.html
# process.ProfilerService = cms.Service (
#         "ProfilerService",
#         firstEvent = cms.untracked.int32(2),
#         lastEvent = cms.untracked.int32(10),
#         paths = cms.untracked.vstring('p1')
# )


# process.GlobalTag.globaltag = "SUPPLY_GLOBAL_TAG"
# process.out.fileName = cms.untracked.string('SUPPLY_OUTPUT_FILE_NAME'),
# process.source.fileNames = cms.untracked.vstring('SUPPLY_INPUT_FILE_NAME')
# process.eventMaker.CMS3tag = cms.string('SUPPLY_CMS3_TAG')
# process.eventMaker.datasetName = cms.string('SUPPLY_DATASETNAME')
# process.maxEvents.input = cms.untracked.int32(SUPPLY_MAX_NEVENTS)

process.GlobalTag.globaltag = "101X_dataRun2_Prompt_v11"
process.out.fileName = cms.untracked.string('ntuple.root')
process.source.fileNames = cms.untracked.vstring('/store/data/Run2018C/MuonEG/MINIAOD/PromptReco-v1/000/319/337/00000/8AA0A1A2-A984-E811-9C77-FA163EA7E2FA.root')
process.eventMaker.CMS3tag = cms.string('SUPPLY_CMS3_TAG')
process.eventMaker.datasetName = cms.string('SUPPLY_DATASETNAME')
process.maxEvents.input = cms.untracked.int32(1000)
