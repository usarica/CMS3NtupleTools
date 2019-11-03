import FWCore.ParameterSet.Config as cms

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#VarParsing_Documentation
# Allow command line options like
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
opts.register('fastsim' , False , mytype=vpbool) # is fastsim?
opts.register('triginfo'  , False , mytype=vpbool) # want (probably broken now) trigger matching information?
opts.register('sparminfo'  , False , mytype=vpbool) # separate flag to enable sparm if fastsim=True isn't specified
opts.register('metrecipe'  , False , mytype=vpbool) # to enable the 2017 94X data,MC MET recipe v2
opts.register('goldenjson'  , "" , mytype=vpstring) # to only process a set of run,lumi sections; see note below for details
opts.register('genxsecanalyzer'  , False , mytype=vpbool) # ONLY run the genxsec analyzer
opts.register('applyEGscalesmear', True , mytype=vpbool) # to enable e/gamma scale and smear corrections
opts.register('dumpAllObjects', False , mytype=vpbool) # if true, use classic edm::Wrapper dumps of the makers
opts.parseArguments()

# this is the section where we try to be a bit smart for the purpose of laziness
if opts.fastsim: opts.data = False
if any(x in opts.inputs for x in ["/store/data/","/Run201"]): opts.data = True
if any(x in opts.inputs for x in ["/store/mc/","MINIAODSIM"]): opts.data = False
if any(x in opts.inputs for x in ["Summer16MiniAODv2","Spring16MiniAODv2","03Feb2017"]): opts.is80x = True
if any(x in opts.inputs for x in ["16MiniAOD","Run2016"]): opts.year = 2016
if any(x in opts.inputs for x in ["17MiniAOD","Run2017"]): opts.year = 2017
if any(x in opts.inputs for x in ["18MiniAOD","Run2018"]): opts.year = 2018

print("""PSet is assuming:
   data? {data} fastsim? {fastsim} is80x? {is80x}
   triginfo? {triginfo} metrecipe? {metrecipe}
   year = {year}
   nevents = {nevents}
   output = {output}
   globaltag = {globaltag}
   inputs = {inputs}
   goldenjson = {goldenjson}
""".format(
   data = opts.data,
   fastsim = opts.fastsim,
   triginfo = opts.triginfo,
   is80x = opts.is80x,
   metrecipe = opts.metrecipe,
   year = opts.year,
   nevents = opts.nevents,
   output = opts.output,
   globaltag = opts.globaltag,
   inputs = opts.inputs,
   goldenjson = opts.goldenjson,
        ))

if not opts.data and opts.year<2016:
   # print(("#"*100+"\n")*2+"[!] This is MC but you've not defined year. To avoid a crash, I'm setting it to 2017\n"+("#"*100+"\n")*2)
   # opts.year=2017
   raise RuntimeError("MC processing must define a year>=2016!")

import CMS3.NtupleMaker.configProcessName as configProcessName
configProcessName.isFastSim=opts.fastsim

# CMS3
process = cms.Process("CMS3")

# Version Control For Python Configuration Files
process.configurationMetadata = cms.untracked.PSet(
        version    = cms.untracked.string('$Revision: 1.11 $'),
        annotation = cms.untracked.string('CMS3'),
        name       = cms.untracked.string('CMS3 test configuration')
)


# load event level configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

# services
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.GlobalTag.globaltag = "94X_mc2017_realistic_v14"
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold  = ''
# hide warning ::getByLabel: An attempt was made to read a Run product before endRun() was called.
process.MessageLogger.suppressWarning = cms.untracked.vstring(["genMaker","sParmMaker"])

process.options = cms.untracked.PSet()


# TODO need to also disable reading this filter decision in the source code for 80x samples
if not opts.is80x:
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
    process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
    baddetEcallist = cms.vuint32(
        [872439604,872422825,872420274,872423218,
         872423215,872416066,872435036,872439336,
         872420273,872436907,872420147,872439731,
         872436657,872420397,872439732,872439339,
         872439603,872422436,872439861,872437051,
         872437052,872420649,872422436,872421950,
         872437185,872422564,872421566,872421695,
         872421955,872421567,872437184,872421951,
         872421694,872437056,872437057,872437313])
    process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
        "EcalBadCalibFilter",
        EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
        ecalMinEt        = cms.double(50.),
        baddetEcal    = baddetEcallist,
        taggingMode = cms.bool(True),
        debug = cms.bool(False)
        )

#Electron Identification for PHYS 14
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupAllVIDIdsInModule, setupVIDElectronSelection, setupVIDPhotonSelection
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
#process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
#process.egmGsfElectronIDSequence = cms.Sequence(process.electronMVAVariableHelper * process.electronMVAValueMapProducer * process.egmGsfElectronIDs)
my_eleid_modules = [
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff'
]
my_phoid_modules = [
    'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff',
    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff'
]

# Load Ntuple producer cff
process.load("CMS3.NtupleMaker.cms3CoreSequences_cff")
if not opts.data: process.load("CMS3.NtupleMaker.cms3GENSequence_cff")
process.load("CMS3.NtupleMaker.cms3PFSequence_cff")
process.eventMaker.isData = cms.bool(opts.data)
process.electronMaker.year = cms.int32(opts.year)
process.photonMaker.year = cms.int32(opts.year)
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

if opts.is80x:
    # DeepCSV computation needed for 80X
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
    deep_discriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfDeepCSVJetTags:probc']
    updateJetCollection( process, btagDiscriminators = deep_discriminators,
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
    )
    updateJetCollection( process, btagDiscriminators = deep_discriminators,
        labelName = 'Puppi',
        jetSource = cms.InputTag('slimmedJetsPuppi'),
        jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
    )
    process.pfJetMaker.pfJetsInputTag = cms.InputTag('selectedUpdatedPatJets')
    process.pfJetPUPPIMaker.pfJetsInputTag = cms.InputTag('selectedUpdatedPatJetsPuppi')
    patAlgosToolsTask = getPatAlgosToolsTask(process)


# Extra trigger information (matching)
if opts.triginfo:
    process.load("CMS3.NtupleMaker.muToTrigAssMaker_cfi")
    process.load("CMS3.NtupleMaker.elToTrigAssMaker_cfi")
    process.muToTrigAssMaker.triggerObjectsName = cms.untracked.string("slimmedPatTrigger")
    process.elToTrigAssMaker.triggerObjectsName = cms.untracked.string("slimmedPatTrigger")
    process.hltMaker.triggerObjectsName = cms.untracked.string("slimmedPatTrigger")
    process.hltMaker.fillTriggerObjects = cms.untracked.bool(True)

if opts.fastsim:
    if opts.is80x:
        process.genMaker.LHEInputTag = cms.InputTag("source")

if opts.year == 2016:
    process.metFilterMaker.doEcalFilterUpdate = False

if opts.is80x and not opts.data:
    process.load("CondCore.CondDB.CondDB_cfi")
    process.CondDB.connect = "frontier://FrontierProd/CMS_CONDITIONS"
    process.l1tPS = cms.ESSource("PoolDBESSource", process.CondDB,
        toGet = cms.VPSet(cms.PSet(record = cms.string("L1TGlobalPrescalesVetosRcd"), tag = cms.string("L1TGlobalPrescalesVetos_passThrough_mc"))))
    process.es_prefer_l1tPS = cms.ESPrefer("PoolDBESSource", "l1tPS")

# Apply E/Gamma corrections if needed
process.load("RecoEgamma.ElectronIdentification.heepIdVarValueMapProducer_cfi")
if (opts.year == 2016):
   if not opts.is80x:
      from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
      setupEgammaPostRecoSeq(process,
                             runEnergyCorrections=opts.applyEGscalesmear,
                             runVID=True,
                             autoAdjustParams=False,
                             eleIDModules = my_eleid_modules,
                             phoIDModules = my_phoid_modules,
                             #Below is from https://github.com/CJLST/ZZAnalysis/blob/Run2Legacy/AnalysisStep/test/MasterPy/ZZ4lAnalysis.py
                             #eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                             #phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                             era='2016-Legacy')
   else:
      # See
      # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#Running_on_2016_MiniAOD_V1
      # and
      # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#Running_on_2017_MiniAOD_V1
      process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
      from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
      setupEgammaPostRecoSeq(process,
                             runEnergyCorrections=opts.applyEGscalesmear,
                             runVID=True,
                             autoAdjustParams=False,
                             eleIDModules = my_eleid_modules,
                             phoIDModules = my_phoid_modules,
                             #Below is from https://github.com/CJLST/ZZAnalysis/blob/Run2Legacy/AnalysisStep/test/MasterPy/ZZ4lAnalysis.py
                             #eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                             #phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                             era='2016-Legacy')


elif (opts.year == 2017):
   from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
   setupEgammaPostRecoSeq(process,
                          runEnergyCorrections=opts.applyEGscalesmear,
                          runVID=True,
                          autoAdjustParams=False,
                          eleIDModules = my_eleid_modules,
                          phoIDModules = my_phoid_modules,
                          #Below is from https://github.com/CJLST/ZZAnalysis/blob/Run2Legacy/AnalysisStep/test/MasterPy/ZZ4lAnalysis.py
                          #phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                          era='2017-Nov17ReReco')

elif (opts.year == 2018):
   from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
   setupEgammaPostRecoSeq(process,
                          runEnergyCorrections=opts.applyEGscalesmear,
                          runVID=True,
                          autoAdjustParams=False,
                          eleIDModules = my_eleid_modules,
                          phoIDModules = my_phoid_modules,
                          #Below is from https://github.com/CJLST/ZZAnalysis/blob/Run2Legacy/AnalysisStep/test/MasterPy/ZZ4lAnalysis.py
                          #eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Autumn18_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                          #phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                          era='2018-Prompt')

#for idmod in my_eleid_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
#for idmod in my_phoid_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

# Setup e/gamma sequence
process.egammaMakerSeq = cms.Sequence( process.heepIDVarValueMaps * process.egammaPostRecoSeq * process.electronMaker * process.photonMaker )

# steal some logic from https://github.com/cms-sw/cmssw/blob/CMSSW_10_4_X/PhysicsTools/NanoAOD/python/nano_cff.py
producers = [
        process.eventMaker,
        process.lumiFilter, # filter after eventmaker so we get run lumi event branches at least, and event counts match up
        process.hltMakerSequence if not opts.fastsim else None,
        process.miniAODrhoSequence,
        process.metFilterMaker,
        process.vertexMaker,
        #process.secondaryVertexMaker,
        process.pfCandidateMaker,
        process.isoTrackMaker,
        process.muonMaker,
        process.egammaMakerSeq,
        process.electronMaker,
        #process.photonMaker,
        process.pfJetMaker,
        process.pfJetPUPPIMaker,
        process.subJetMaker,
        process.pfmetMaker,
        process.pfmetpuppiMaker,
        process.pftauMaker,
        # Disable these temporarily until they are renewed
        #process.muToTrigAssMaker if opts.triginfo else None,
        #process.elToTrigAssMaker if opts.triginfo else None,
        process.genMaker if not opts.data else None,
        process.genJetMaker if not opts.data else None,
        #process.candToGenAssMaker if not opts.data else None,
        #process.pdfinfoMaker if not opts.data else None,
        process.puSummaryInfoMaker if not opts.data else None,
        #process.hypDilepMaker,
        #process.sParmMaker if (opts.fastsim or opts.sparminfo) else None,
        ]
if opts.genxsecanalyzer and not opts.data:
    process.genxsecanalyzer = cms.EDAnalyzer("GenXSecAnalyzer")
    producers = [process.genxsecanalyzer]
total_path = None
for ip,producer in enumerate(producers):
    if producer is None: continue
    if ip == 0:
        total_path = producer
        continue

    if opts.is80x and producer in [process.isoTrackMaker]: continue

    if opts.fastsim:
        if opts.is80x and producer in [process.metFilterMaker]: continue

    if not opts.is80x and producer in [process.metFilterMaker]: total_path *= process.ecalBadCalibReducedMINIAODFilter

    if opts.metrecipe and producer == process.pfmetMaker:
        total_path *= process.fullPatMetSequenceModifiedMET * process.pfmetMaker * process.pfmetMakerModifiedMET
        continue

    if producer == process.electronMaker:
    #   total_path *= process.heepIDVarValueMaps * process.egammaPostRecoSeq * process.electronMaker
       total_path *= process.egammaMakerSeq
       continue


    total_path *= producer

process.p = cms.Path(total_path)

if opts.is80x:
    from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
    process.p = cms.Path(process.p._seq,patAlgosToolsTask)

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

# Event maker
process.eventMaker.CMS3tag = cms.string('SUPPLY_CMS3_TAG')
process.eventMaker.datasetName = cms.string('SUPPLY_DATASETNAME')

# Final output
if opts.dumpAllObjects:
   process.out = cms.OutputModule(
      "PoolOutputModule",
      fileName     = cms.untracked.string(opts.output),
      dropMetaData = cms.untracked.string("ALL"),
      basketSize = cms.untracked.int32(16384*23)
   )
   process.out.outputCommands = cms.untracked.vstring( 'drop *' )
   process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS3*'))
   process.outpath = cms.EndPath(process.out)

else:
   process.TFileService = cms.Service(
      'TFileService',
      fileName=cms.string(opts.output)
   )

   process.load("CMS3.NtupleMaker.cms3Ntuplizer_cfi")
   process.cms3ntuple.year = cms.int32(opts.year)
   process.cms3ntuple.isMC = cms.bool((not opts.data))
   process.outpath = cms.EndPath(process.cms3ntuple)


fprocdump = open(opts.output.replace('.root','_run_cfg.py'),'w')
fprocdump.write(process.dumpPython())
fprocdump.close()
