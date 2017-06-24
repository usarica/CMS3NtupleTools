import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff   import *


is_data = False
is_prompt = False
is_fastsim = False
is_relval = False

do_deepbtag = True

import CMS3.NtupleMaker.configProcessName as configProcessName
configProcessName.name="PAT"
if is_data and is_prompt:
    configProcessName.name="RECO"

configProcessName.name2="RECO"

if is_relval:
    configProcessName.name="reRECO"
    configProcessName.name2="reRECO"

if is_fastsim:
    configProcessName.fastSimName="HLT"
    configProcessName.name2=configProcessName.fastSimName
configProcessName.isFastSim=is_fastsim

# CMS3
process = cms.Process("CMS3")

# Version Control For Python Configuration Files
process.configurationMetadata = cms.untracked.PSet(
        version    = cms.untracked.string('$Revision: 1.11 $'),
        annotation = cms.untracked.string('CMS3'),
        name       = cms.untracked.string('CMS3 test configuration')
)

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
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True),SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName     = cms.untracked.string('ntuple.root'),
                               # fileName     = cms.untracked.string('ntuple2.root'),
                               dropMetaData = cms.untracked.string("ALL"),
                               # basketSize = cms.untracked.int32(16384*150)
                               # basketSize = cms.untracked.int32(16384*10)
                               basketSize = cms.untracked.int32(16384*23)
)



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
                 ]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# Load Ntuple producer cff
process.load("CMS3.NtupleMaker.cms3CoreSequences_cff")
if not is_data: process.load("CMS3.NtupleMaker.cms3GENSequence_cff")
process.load("CMS3.NtupleMaker.cms3PFSequence_cff")
process.eventMaker.isData                        = cms.bool(is_data)
    
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
applyResiduals=is_data #application of residual corrections. Have to be set to True once the 13 TeV residual corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or for analysis checks and developments (not the official recommendation!).
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

if not is_data:
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    #default configuration for miniAOD reprocessing, change the isData flag to run on data
    #for a full met computation, remove the pfCandColl input
    runMetCorAndUncFromMiniAOD(process,
                               isData=is_data,
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

if is_data:
    process.p = cms.Path( 
        process.metFilterMaker *
        process.egmGsfElectronIDSequence *     
        process.vertexMaker *
        process.secondaryVertexMaker *
        process.eventMaker *
        process.pfCandidateMaker *
        #  process.isoTrackMaker *
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
        #  process.isoTrackMaker *
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

# process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v7"
# process.out.fileName = cms.untracked.string('ntuple.root')
# process.source.fileNames = cms.untracked.vstring('file:DataDoubleEG2016C.root')
# process.eventMaker.CMS3tag = cms.string('V08-00-18')
# process.eventMaker.datasetName = cms.string('/DoubleEG/Run2016C-03Feb2017-v1/MINIAOD')
# process.maxEvents.input = cms.untracked.int32(3000)
