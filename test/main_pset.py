import FWCore.ParameterSet.Config as cms

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#VarParsing_Documentation
# Allow command line options like
#     cmsRun main_pset.py data=True prompt=True   # prompt data
#     cmsRun main_pset.py data=False               # MC
#     cmsRun main_pset.py fastsim=True             # fastsim
import FWCore.ParameterSet.VarParsing as VarParsing
opts = VarParsing.VarParsing('python')
vpbool = VarParsing.VarParsing.varType.bool
vpint = VarParsing.VarParsing.varType.int
vpstring = VarParsing.VarParsing.varType.string
opts.register('data'    , False  , mytype=vpbool)
opts.register('globaltag'    , ""  , mytype=vpstring)
opts.register('inputs'    , ""  , mytype=vpstring) # comma separated list of input files
opts.register('output'    , "ntuple.root"  , mytype=vpstring)
opts.register('nevents'    , -1  , mytype=vpint)
opts.register('year'    , -1  , mytype=vpint) # year for MC weight and other purposes (2016,2017,2018); defaults to 2017
opts.register('is80x'    , False  , mytype=vpbool) # is 2016 80X sample?
opts.register('prompt'  , False  , mytype=vpbool)
opts.register('fastsim' , False , mytype=vpbool)
opts.register('relval'  , False , mytype=vpbool)
opts.register('triginfo'  , False , mytype=vpbool)
opts.register('metrecipe'  , False , mytype=vpbool) # to enable the 2017 94X data,MC MET recipe v2
opts.register('eventmakeronly'  , False , mytype=vpbool)
opts.register('goldenjson'  , "" , mytype=vpstring) # to only process a set of run,lumi sections; see note below for details
opts.register('name'  , "" , mytype=vpstring) # hacky variable to override name for samples where last path/process is "DQM"
opts.parseArguments()
# be smart. if fastsim, it's obviously MC
# if it's MC, it's obviously not prompt
if opts.fastsim: opts.data = False
if not opts.data: opts.prompt = False
print("""PSet is assuming:
   data? {data} prompt? {prompt} fastsim? {fastsim} relval? {relval} is80x? {is80x}
   triginfo? {triginfo} metrecipe? {metrecipe} eventmakeronly? {eventmakeronly}
   year = {year}
   nevents = {nevents}
   output = {output}
   name = {name}
   globaltag = {globaltag}
   inputs = {inputs}
   goldenjson = {goldenjson}
""".format(
   data = opts.data,
   prompt = opts.prompt,
   fastsim = opts.fastsim,
   relval = opts.relval,
   triginfo = opts.triginfo,
   is80x = opts.is80x,
   metrecipe = opts.metrecipe,
   eventmakeronly = opts.eventmakeronly,
   year = opts.year,
   nevents = opts.nevents,
   output = opts.output,
   name = opts.name,
   globaltag = opts.globaltag,
   inputs = opts.inputs,
   goldenjson = opts.goldenjson,
        ))

if not opts.data and opts.year<2016:
    print(("#"*100+"\n")*2+"[!] This is MC but you've not defined year. To avoid a crash, I'm setting it to 2017\n"+("#"*100+"\n")*2)
    opts.year=2017
  # raise RuntimeError("MC processing must define a year>=2016!")

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

if str(opts.name).strip():
    configProcessName.name = str(opts.name).strip()

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
process.GlobalTag.globaltag = "94X_mc2017_realistic_v14"
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold  = ''
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(False) )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName     = cms.untracked.string('ntuple.root'),
                               dropMetaData = cms.untracked.string("ALL"),
                               basketSize = cms.untracked.int32(16384*23)
)

print "FIXME do we need this?"
print "FIXME do we need this?"
print "FIXME do we need this?"
print "FIXME do we need this?"
print "FIXME do we need this?"
print "FIXME do we need this?"
print "FIXME do we need this?"
print "FIXME do we need this?"
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
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

#Electron Identification for PHYS 14
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupAllVIDIdsInModule, setupVIDElectronSelection
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons',"",configProcessName.name)
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons',"",configProcessName.name)
process.egmGsfElectronIDSequence = cms.Sequence(process.electronMVAVariableHelper * process.electronMVAValueMapProducer * process.egmGsfElectronIDs)
my_id_modules = [
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
                 ]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# Load Ntuple producer cff
process.load("CMS3.NtupleMaker.cms3CoreSequences_cff")
if not opts.data: process.load("CMS3.NtupleMaker.cms3GENSequence_cff")
process.load("CMS3.NtupleMaker.cms3PFSequence_cff")
process.eventMaker.isData                        = cms.bool(opts.data)
if not opts.data:
    process.genMaker.year = cms.int32(opts.year)

#Options for Input
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'file:C6BB52E8-F341-E811-8A2F-001E677927EC.root',
                            )
)
def find_up(fname):
    import os
    d = os.getcwd()
    while d != "/":
        t, d = os.path.join(d,fname), os.path.dirname(d)
        if os.path.exists(t): return t

if opts.goldenjson and find_up(opts.goldenjson):
    goldenjson = find_up(opts.goldenjson)
    # if we filter in the process.source, then the events are just skipped
    # so we use a custom lumiFilter to skip *after* the EventMaker to keep
    # total event counts in agreement with DBS, but also have evt_event,run,lumiBlock
    # for babymakers to filter
    skip_event = False
    import FWCore.PythonUtilities.LumiList as LumiList
    # JSONfile = "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"
    lumilist = LumiList.LumiList(filename=goldenjson).getCMSSWString().split(',')
    print("Found json list of lumis to process with {} lumi sections from {}".format(len(lumilist),goldenjson))
    print("Skipping {} if they're not in the lumi list".format("events entirely" if skip_event else "anything after eventMaker"))
    if skip_event:
        process.source.lumisToProcess = cms.untracked(cms.VLuminosityBlockRange()+lumilist)
    else:
        process.lumiFilter.lumisToProcess = cms.untracked(cms.VLuminosityBlockRange()+lumilist)

#Max Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(opts.nevents) )

process.outpath = cms.EndPath(process.out)
process.out.outputCommands = cms.untracked.vstring( 'drop *' )

extra = {}
if opts.metrecipe:
    process.pfmetMakerModifiedMET = process.pfmetMaker.clone()
    process.pfmetMakerModifiedMET.pfMetInputTag_ = cms.InputTag("slimmedMETsModifiedMET","","CMS3")
    process.pfmetMaker.aliasPrefix = cms.untracked.string("evt_old")
    extra = dict(
            fixEE2017=True,
            fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139},
            postfix="ModifiedMET",
            )

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=opts.data,
                           **extra
                           )

do_deepbtag = opts.is80x
if do_deepbtag:
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    deep_discriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfDeepCSVJetTags:probc']
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
       jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
        btagDiscriminators = deep_discriminators
    )
    updateJetCollection(
        process,
        labelName = 'Puppi',
        jetSource = cms.InputTag('slimmedJetsPuppi'),
       jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
        btagDiscriminators = deep_discriminators
    )
    process.pfJetMaker.pfJetsInputTag = cms.InputTag('selectedUpdatedPatJets')
    process.pfJetPUPPIMaker.pfJetsInputTag = cms.InputTag('selectedUpdatedPatJetsPuppi')

process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS3*'))

### -------------------------------------------------------------------
### the lines below remove the L2L3 residual corrections when processing MC
### -------------------------------------------------------------------
#Run corrected MET maker
applyResiduals=opts.data #application of residual corrections. Have to be set to True once the 13 TeV residual corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or for analysis checks and developments (not the official recommendation!).
if not applyResiduals:
    process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    # process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    # process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
### ------------------------------------------------------------------

# end Run corrected MET maker

# Extra trigger information (matching)
if opts.triginfo:
    process.load("CMS3.NtupleMaker.muToTrigAssMaker_cfi")
    process.load("CMS3.NtupleMaker.elToTrigAssMaker_cfi")
    if opts.data and opts.prompt:
        process.muToTrigAssMaker.triggerObjectsName = cms.untracked.string("slimmedPatTrigger")
        process.elToTrigAssMaker.triggerObjectsName = cms.untracked.string("slimmedPatTrigger")
        process.hltMaker.triggerObjectsName = cms.untracked.string("slimmedPatTrigger")
    process.hltMaker.fillTriggerObjects = cms.untracked.bool(True)

# steal some logic from https://github.com/cms-sw/cmssw/blob/CMSSW_10_4_X/PhysicsTools/NanoAOD/python/nano_cff.py
producers = [
        process.eventMaker,
        process.lumiFilter, # filter after eventmaker so we get run lumi event branches at least, and event counts match up
        process.metFilterMaker,
        process.egmGsfElectronIDSequence,
        process.vertexMaker,
        process.secondaryVertexMaker,
        process.pfCandidateMaker,
        process.isoTrackMaker,
        process.electronMaker,
        process.muonMaker,
        process.pfJetMaker,
        process.pfJetPUPPIMaker,
        process.subJetMaker,
        process.pfmetMaker,
        process.pfmetpuppiMaker,
        process.hltMakerSequence,
        process.miniAODrhoSequence,
        process.pftauMaker,
        process.photonMaker,
        process.muToTrigAssMaker if opts.triginfo else None,
        process.elToTrigAssMaker if opts.triginfo else None,
        process.genMaker if not opts.data else None,
        process.genJetMaker if not opts.data else None,
        process.candToGenAssMaker if not opts.data else None,
        process.pdfinfoMaker if not opts.data else None,
        process.puSummaryInfoMaker if not opts.data else None,
        process.hypDilepMaker,
        ]
total_path = None
for ip,producer in enumerate(producers):
    if producer is None: continue
    if ip == 0:
        total_path = producer
        continue

    if opts.is80x and producer in [process.isoTrackMaker]:
        continue

    if opts.eventmakeronly and producer not in [
            process.eventMaker
            ]: continue

    if opts.metrecipe and producer == process.pfmetMaker:
        total_path *= process.fullPatMetSequenceModifiedMET * process.pfmetMaker * process.pfmetMakerModifiedMET
        continue

    total_path *= producer

process.p = cms.Path(total_path)

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

if not opts.globaltag:
    process.GlobalTag.globaltag = "102X_upgrade2018_realistic_v15"
else:
    process.GlobalTag.globaltag = opts.globaltag
if not opts.inputs:
    process.source.fileNames = cms.untracked.vstring("/store/mc/RunIIAutumn18MiniAOD/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/84347EC0-60B4-5145-8F92-37F1975CA79D.root")
else:
    inputs = opts.inputs.split(",")
    print("Running on inputs: {}".format(inputs))
    process.source.fileNames = cms.untracked.vstring(inputs)
process.out.fileName = cms.untracked.string(opts.output)
process.eventMaker.CMS3tag = cms.string('SUPPLY_CMS3_TAG')
process.eventMaker.datasetName = cms.string('SUPPLY_DATASETNAME')
