import FWCore.ParameterSet.Config as cms

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#VarParsing_Documentation
# Allow command line options like
#     cmsRun main_pset.py data=True prompt=True   # prompt data
#     cmsRun main_pset.py data=False               # MC
#     cmsRun main_pset.py fastsim=True             # fastsim
import FWCore.ParameterSet.VarParsing as VarParsing
opts = VarParsing.VarParsing('python')
vpbool = VarParsing.VarParsing.varType.bool
opts.register('data'    , True  , mytype=vpbool)
opts.register('prompt'  , True  , mytype=vpbool)
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
process.GlobalTag.globaltag = "80X_mcRun2_asymptotic_2016_miniAODv2_v0" #80X
#process.GlobalTag.globaltag = "91X_upgrade2017_realistic_v5" #MC
#process.GlobalTag.globaltag = "91X_dataRun2_relval_v6" #data
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
process.egmGsfElectronIDSequence = cms.Sequence(process.electronMVAValueMapProducer * process.egmGsfElectronIDs)
my_id_modules = [
        # 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
        # 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
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
                                # 'file:/hadoop/cms/phedex/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/7AEAFCAD-266F-E511-8A2A-001E67A3F3DF.root',
                                # 'root://cmsxrootd.fnal.gov//store/mc/RunIIFall15MiniAODv1/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/0E47EC63-7B9D-E511-B714-B083FED426E5.root
#         'file:/hadoop/cms/phedex/store/mc/RunIISpring16MiniAODv1/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/60000/F4EA8D09-9002-E611-9D1B-1CC1DE19274E.root',
#                                '/store/mc/RunIISpring16MiniAODv2/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/D63C4E53-D91B-E611-AC83-FA163E5810F7.root',
                                # 'file:RelValProdQCD_Pt_3000_3500_13.root'
                                # 'file:/home/users/namin/2017/slimming/CMSSW_8_0_26_patch1/src/CMS3/NtupleMaker/test/A8B84A69-C1D7-E611-831F-5065F382B2D1.root',
                                # 'file:/home/users/namin/2017/slimming/CMSSW_8_0_26_patch1/src/CMS3/NtupleMaker/test/A8B84A69-C1D7-E611-831F-5065F382B2D1.root',
                                # 'file:/home/users/namin/2017/slimming/CMSSW_8_0_26_patch1/src/CMS3/NtupleMaker/test',
                                #'file:/home/users/namin/2017/slimming/CMSSW_8_0_26_patch1/src/CMS3/NtupleMaker/test/TTJets_HT-1200to2500.root',
                                # 'file:DataDoubleEG2016C.root',
                                # 'file:QCD_HT200to300.root',
                                # 'file:20457CC1-74D7-E611-A445-24BE05CE2E81.root',

#                                'root://cmsxrootd.fnal.gov//store/relval/CMSSW_9_2_0/SingleMuon/MINIAOD/91X_dataRun2_relval_v6_RelVal_sigMu2016B-v1/10000/746430BE-773C-E711-8419-0CC47A745298.root',
                                #'root://cmsxrootd.fnal.gov//store/relval/CMSSW_9_2_0/SingleMuon/MINIAOD/91X_dataRun2_relval_v6_RelVal_sigMu2016E-v1/10000/5C79F5F3-B13C-E711-AEFD-0CC47A4D762A.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/PhaseISpring17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/FlatPU28to62_90X_upgrade2017_realistic_v20-v1/00000/02781287-E22A-E711-8EF8-A0000420FE80.root',
                                'file:/home/users/mderdzinski/ntupling/CMSSW_8_0_26_patch1_CMS4_V00-00-02/src/CMS3/NtupleMaker/TTJets_HT-1200to2500.root',
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

# process.GlobalTag.globaltag = "94X_dataRun2_ReReco_EOY17_v2"
process.out.fileName = cms.untracked.string('ntuple.root')
# process.source.fileNames = cms.untracked.vstring('file:/home/users/namin/2017/lepmvacms4/CMSSW_9_4_0/src/CMS3/NtupleMaker/test/EAED912B-F7DE-E711-8E9B-0242AC1C0500.root')
# process.source.fileNames = cms.untracked.vstring('/store/data/Run2017F/DoubleEG/MINIAOD/17Nov2017-v1/60000/EAED912B-F7DE-E711-8E9B-0242AC1C0500.root')
# process.source.fileNames = cms.untracked.vstring('/store/data/Run2017E/HTMHT/MINIAOD/31Mar2018-v1/90000/D2735DEC-0D37-E811-AE40-A4BF0115947C.root')
process.source.fileNames = cms.untracked.vstring('/store/data/Run2017D/SingleMuon/MINIAOD/31Mar2018-v1/80000/1E703527-F436-E811-80A7-E0DB55FC1055.root')
# process.source.fileNames = cms.untracked.vstring('file:1E703527-F436-E811-80A7-E0DB55FC1055.root')
process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v6'
process.eventMaker.CMS3tag = cms.string('blah')
process.eventMaker.datasetName = cms.string('blah')
process.maxEvents.input = cms.untracked.int32(1000)
