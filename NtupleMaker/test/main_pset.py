import FWCore.ParameterSet.Config as cms

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#VarParsing_Documentation
# Allow command line options like
#     cmsRun main_pset.py data=False               # MC
#     cmsRun main_pset.py fastsim=True             # fastsim
import FWCore.ParameterSet.VarParsing as VarParsing
opts = VarParsing.VarParsing('python')
vpbool = VarParsing.VarParsing.varType.bool
vpint = VarParsing.VarParsing.varType.int
vpfloat = VarParsing.VarParsing.varType.float
vpstring = VarParsing.VarParsing.varType.string
opts.register('dumpProcess'    , False  , mytype=vpbool)
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
opts.register('applyMuoncorr', True , mytype=vpbool) # to enable muon scale and smear corrections
opts.register('updatePileupJetId', True , mytype=vpbool) # to enable dating the pile-up jet id
opts.register('keepGenParticles' , "reducedfinalstates" , mytype=vpstring) # to keep gen. particles. See CMS3NtupleMaker::ParticleRecordLevel enums
opts.register('keepGenJets' , True , mytype=vpbool) # to keep gen. jets
opts.register('dumpAllObjects', False , mytype=vpbool) # if true, use classic edm::Wrapper dumps of the makers
opts.register('xsec', -1, mytype=vpfloat) # xsec value of the MC sample in pb, hopefully
opts.register('BR', -1, mytype=vpfloat) # BR value of the MC sample
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

# Load the conditions database
process.load("CondCore.CondDB.CondDB_cfi")

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
# Add 2016 and 2018 HZZ MVA retraining
if opts.year == 2016:
   my_eleid_modules.append("RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff")
# 2017 uses Fall17V2_Iso
elif opts.year == 2018:
   my_eleid_modules.append("RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Autumn18_ID_ISO_cff")
my_phoid_modules = [
   'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff',
   'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff'
   ]

# Load Ntuple producer cff
process.load("CMS3.NtupleMaker.cms3CoreSequences_cff")
if not opts.data: process.load("CMS3.NtupleMaker.cms3GENSequence_cff")
#process.eventMaker.isData = cms.bool(opts.data)
process.muonMaker.year = cms.int32(opts.year)
process.muonMaker.refurbishSelections = cms.bool(opts.is80x and not opts.data)
process.electronMaker.year = cms.int32(opts.year)
process.photonMaker.year = cms.int32(opts.year)
process.isoTrackMaker.year = cms.int32(opts.year)
if not opts.data:
   process.genMaker.year = cms.int32(opts.year)
   process.genMaker.xsec = cms.double(opts.xsec)
   process.genMaker.BR = cms.double(opts.BR)

#Options for Input
process.source = cms.Source(
   "PoolSource",
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
   return None

if opts.data and opts.goldenjson:
    goldenjson = find_up(opts.goldenjson)
    if goldenjson is None:
       goldenjson = find_up('data/LumiJSON/'+opts.goldenjson)
    if goldenjson is None:
       raise RuntimeError("Golden JSON file {} cannot be found!".format(opts.goldenjson))
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

# Need to load L1 prescales in 8.0.X MC
if opts.is80x and not opts.data:
   process.CondDB.connect = "frontier://FrontierProd/CMS_CONDITIONS"
   process.l1tPS = cms.ESSource(
      "PoolDBESSource", process.CondDB,
      toGet = cms.VPSet(
         cms.PSet(
            record = cms.string("L1TGlobalPrescalesVetosRcd"),
            tag = cms.string("L1TGlobalPrescalesVetos_passThrough_mc")
            )
         )
      )
   process.es_prefer_l1tPS = cms.ESPrefer("PoolDBESSource", "l1tPS")

# L1 prefiring issue for 2016 and 2017 data
## Recipe taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe#Recipe_details_10_2_X_X_10_or_9
prefiringWeightsTag=""
if not opts.data and (opts.year == 2016 or opts.year == 2017):
   from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
   prefiringdataera=""
   if opts.year == 2016:
      prefiringdataera="2016BtoH"
   else:
      prefiringdataera="2017BtoF"
   prefiringWeightsTag = "prefiringweight"
   process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
      DataEra = cms.string(prefiringdataera),
      UseJetEMPt = cms.bool(False),
      PrefiringRateSystematicUncty = cms.double(0.2),
      SkipWarnings = False
      )
   process.Prefiring = cms.Path(process.prefiringweight)


# Update jet variables
## Pile-up jet id
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
    jets=cms.InputTag("slimmedJets"),
    inputIsCorrected=False,
    applyJec=True,
    vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
    )
ak4jetsTag="AK4PFchs"
ak4puppijetsTag="AK4PFPuppi"
ak8jetsTag="AK8PFPuppi"
if opts.is80x:
   ak8jetsTag="AK8PFchs"
## Update the jet collection tags
process.pfJetMaker.jetCollection = cms.untracked.string(ak4jetsTag)
process.pfJetPUPPIMaker.jetCollection = cms.untracked.string(ak4puppijetsTag)
process.subJetMaker.jetCollection = cms.untracked.string(ak8jetsTag)
process.pfJetMaker.isMC = cms.bool((not opts.data))
process.pfJetPUPPIMaker.isMC = cms.bool((not opts.data))
process.subJetMaker.isMC = cms.bool((not opts.data))
## Reapply JECs
### Taken from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
jecVersion=""
if opts.year == 2016:
   if opts.data:
      jecVersion = "Summer16_07Aug2017All_V11_DATA"
   else:
      jecVersion = "Summer16_07Aug2017_V11_MC"
elif opts.year == 2017:
   if opts.data:
      jecVersion = "Fall17_17Nov2017_V32_102X_DATA"
   else:
      jecVersion = "Fall17_17Nov2017_V32_102X_MC"
elif opts.year == 2018:
   if opts.data:
      jecVersion = "Autumn18_RunABCD_V19_DATA"
   else:
      jecVersion = "Autumn18_V19_MC"
if jecVersion != "":
   process.jec = cms.ESSource(
      "PoolDBESSource",
      DBParameters = cms.PSet(
         messageLevel = cms.untracked.int32(1)
         ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
         cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_{}_{}'.format(jecVersion, ak4jetsTag)),
            label  = cms.untracked.string(ak4jetsTag)
            ),
         cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_{}_{}'.format(jecVersion, ak4puppijetsTag)),
            label  = cms.untracked.string(ak4puppijetsTag)
            ),
         cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_{}_{}'.format(jecVersion, ak8jetsTag)),
            label  = cms.untracked.string(ak8jetsTag)
            ),
         ),
      connect = cms.string('sqlite_fip:CMS3/NtupleMaker/data/JECs/{}.db'.format(jecVersion)),
      )

   ## Add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
   process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

   from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
   from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets

   ###################
   ### AK4 PF JETS ###
   ###################
   ## Reapply JECs here
   process.slimmedJetCorrFactors = updatedPatJetCorrFactors.clone(
      src = cms.InputTag("slimmedJets"),
      primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
      levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute']),
      payload = ak4jetsTag
      )
   ## Data applies L2L3Residual corrections as well
   if opts.data:
      process.slimmedJetCorrFactors.levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])
   ## This is the new input jet collection
   process.slimmedCorrectedJets = updatedPatJets.clone(
      jetSource = cms.InputTag("slimmedJets"),
      jetCorrFactorsSource = cms.VInputTag(cms.InputTag("slimmedJetCorrFactors"))
      )
   ## Add pileup id and discriminant to patJetsReapplyJEC
   process.slimmedCorrectedJets.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
   process.slimmedCorrectedJets.userData.userInts.src += ['pileupJetIdUpdated:fullId']
   ## Replace inputs from slimmedJets
   #process.QGTagger.srcJets = cms.InputTag('slimmedCorrectedJets')
   process.pfJetMaker.pfJetsInputTag = cms.InputTag('slimmedCorrectedJets')

   ######################
   ### AK4 PUPPI JETS ###
   ######################
   ## Reapply JECs here
   process.slimmedPuppiJetCorrFactors = updatedPatJetCorrFactors.clone(
      src = cms.InputTag("slimmedJetsPuppi"),
      primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
      levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute']),
      payload = ak4puppijetsTag
      )
   ## Data applies L2L3Residual corrections as well
   if opts.data:
      process.slimmedPuppiJetCorrFactors.levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])
   ## This is the new input jet collection
   process.slimmedCorrectedJetsPuppi = updatedPatJets.clone(
      jetSource = cms.InputTag("slimmedJetsPuppi"),
      jetCorrFactorsSource = cms.VInputTag(cms.InputTag("slimmedPuppiJetCorrFactors"))
      )
   ## Replace inputs from slimmedJets
   process.pfJetPUPPIMaker.pfJetsInputTag = cms.InputTag('slimmedCorrectedJetsPuppi')

   ################
   ### AK8 JETS ###
   ################
   ## Reapply JECs here
   process.slimmedJetAK8CorrFactors = updatedPatJetCorrFactors.clone(
      src = cms.InputTag("slimmedJetsAK8"),
      primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
      levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute']),
      payload = ak8jetsTag
      )
   ## Data applies L2L3Residual corrections as well
   if opts.data:
      process.slimmedJetAK8CorrFactors.levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])
   ## This is the new input jet collection
   process.slimmedCorrectedJetsAK8 = updatedPatJets.clone(
      jetSource = cms.InputTag("slimmedJetsAK8"),
      jetCorrFactorsSource = cms.VInputTag(cms.InputTag("slimmedJetAK8CorrFactors"))
      )
   ## Replace inputs from slimmedJets
   process.subJetMaker.pfJetsInputTag = cms.InputTag('slimmedCorrectedJetsAK8')
else:
   raise RuntimeError("JEC version is unknown!")
## Re-apply JERs
### Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
jerVersion=""
if opts.year == 2016:
   if opts.data:
      jerVersion = "Summer16_25nsV1_DATA"
   else:
      jerVersion = "Summer16_25nsV1_MC"
elif opts.year == 2017:
   if opts.data:
      jerVersion = "Fall17_V3_102X_DATA"
   else:
      jerVersion = "Fall17_V3_102X_MC"
elif opts.year == 2018:
   if opts.data:
      jerVersion = "Autumn18_V7b_DATA"
   else:
      jerVersion = "Autumn18_V7b_MC"
if jerVersion != "":
   process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
   process.jer = cms.ESSource(
      "PoolDBESSource",
      DBParameters = cms.PSet(
         messageLevel = cms.untracked.int32(1)
         ),
      #timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
         cms.PSet(
            record = cms.string('JetResolutionRcd'),
            tag    = cms.string('JR_{}_PtResolution_{}'.format(jerVersion, ak4jetsTag)),
            label  = cms.untracked.string("{}_pt".format(ak4jetsTag))
            ),
         #cms.PSet(
            #record = cms.string('JetResolutionRcd'),
            #tag    = cms.string('JR_{}_EtaResolution_{}'.format(jerVersion, ak4jetsTag)),
            #label  = cms.untracked.string("{}_eta".format(ak4jetsTag))
            #),
         #cms.PSet(
            #record = cms.string('JetResolutionRcd'),
            #tag    = cms.string('JR_{}_PhiResolution_{}'.format(jerVersion, ak4jetsTag)),
            #label  = cms.untracked.string("{}_phi".format(ak4jetsTag))
            #),
         cms.PSet(
            record = cms.string('JetResolutionScaleFactorRcd'),
            tag    = cms.string('JR_{}_SF_{}'.format(jerVersion, ak4jetsTag)),
            label  = cms.untracked.string(ak4jetsTag)
            ),
         cms.PSet(
            record = cms.string('JetResolutionRcd'),
            tag    = cms.string('JR_{}_PtResolution_{}'.format(jerVersion, ak4puppijetsTag)),
            label  = cms.untracked.string("{}_pt".format(ak4puppijetsTag))
            ),
         #cms.PSet(
            #record = cms.string('JetResolutionRcd'),
            #tag    = cms.string('JR_{}_EtaResolution_{}'.format(jerVersion, ak4puppijetsTag)),
            #label  = cms.untracked.string("{}_eta".format(ak4puppijetsTag))
            #),
         #cms.PSet(
            #record = cms.string('JetResolutionRcd'),
            #tag    = cms.string('JR_{}_PhiResolution_{}'.format(jerVersion, ak4puppijetsTag)),
            #label  = cms.untracked.string("{}_phi".format(ak4puppijetsTag))
            #),
         cms.PSet(
            record = cms.string('JetResolutionScaleFactorRcd'),
            tag    = cms.string('JR_{}_SF_{}'.format(jerVersion, ak4puppijetsTag)),
            label  = cms.untracked.string(ak4puppijetsTag)
            ),
         cms.PSet(
            record = cms.string('JetResolutionRcd'),
            tag    = cms.string('JR_{}_PtResolution_{}'.format(jerVersion, ak8jetsTag)),
            label  = cms.untracked.string("{}_pt".format(ak8jetsTag))
            ),
         #cms.PSet(
            #record = cms.string('JetResolutionRcd'),
            #tag    = cms.string('JR_{}_EtaResolution_{}'.format(jerVersion, ak8jetsTag)),
            #label  = cms.untracked.string("{}_eta".format(ak8jetsTag))
            #),
         #cms.PSet(
            #record = cms.string('JetResolutionRcd'),
            #tag    = cms.string('JR_{}_PhiResolution_{}'.format(jerVersion, ak8jetsTag)),
            #label  = cms.untracked.string("{}_phi".format(ak8jetsTag))
            #),
         cms.PSet(
            record = cms.string('JetResolutionScaleFactorRcd'),
            tag    = cms.string('JR_{}_SF_{}'.format(jerVersion, ak8jetsTag)),
            label  = cms.untracked.string(ak8jetsTag)
            ),
         ),
      connect = cms.string('sqlite_fip:CMS3/NtupleMaker/data/JERs/{}.db'.format(jerVersion)),
      )
   process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')
elif not opts.data:
   raise RuntimeError("JER version is unknown!")



# MET corrections and uncertainties
## Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_or_10
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
## PF MET
extra = dict(
   postfix = "ModifiedMET"
   )
### MET v2 recipe for 2017 PF MET
if opts.metrecipe:
   extra = dict(
      fixEE2017 = True,
      fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139},
      postfix = "ModifiedMET",
      )
process.pfmetMaker.metSrc = cms.InputTag("slimmedMETsModifiedMET","","CMS3")
runMetCorAndUncFromMiniAOD(
   process,
   isData=opts.data,
   **extra
   )
## PUPPI MET
makePuppiesFromMiniAOD( process, True ) # This function resets photon IDs!!!
extra_puppi = dict(
   metType = "Puppi",
   jetFlavor = ak4puppijetsTag,
   postfix = "ModifiedPuppiMET",
   )
### MET v2 recipe for 2017 here as well, but it is not clear whether this is official or not
if opts.metrecipe:
   extra_puppi = dict(
      metType = "Puppi",
      jetFlavor = ak4puppijetsTag,
      fixEE2017 = True,
      fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139},
      postfix = "ModifiedPuppiMET",
      )
process.pfmetpuppiMaker.metSrc = cms.InputTag("slimmedMETsModifiedPuppiMET","","CMS3")
runMetCorAndUncFromMiniAOD(
   process,
   isData=opts.data,
   **extra_puppi
   )
### Have no idea about what these things do, but safer to recompute weights I suppose...
process.puppiNoLep.useExistingWeights = False
process.puppi.useExistingWeights = False
## These variables are somehow dropped
#process.slimmedCorrectedJets.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
#process.slimmedCorrectedJets.userData.userInts.src += ['pileupJetIdUpdated:fullId']


# Apply E/Gamma corrections if needed
## Do this after re-applying JECs on Puppi MET!!!
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

from CMS3.NtupleMaker.utils.replaceMVAValuesByRaw import getMVACutflowDictionary,replaceMVAValuesByRaw # Needed to get Rawvalues instead of values
allMVAcuts=cms.VPSet()
electronMVAcuts=cms.VPSet()
photonMVAcuts=cms.VPSet()
if hasattr(process,"egmGsfElectronIDs"):
   electronMVAcuts = getMVACutflowDictionary(process.egmGsfElectronIDs)
if hasattr(process,"egmPhotonIDs"):
   photonMVAcuts = getMVACutflowDictionary(process.egmPhotonIDs)
for p in electronMVAcuts:
   allMVAcuts.append(p)
for p in photonMVAcuts:
   allMVAcuts.append(p)
replaceMVAValuesByRaw(process.slimmedElectrons, allMVAcuts)
replaceMVAValuesByRaw(process.slimmedPhotons, allMVAcuts)
process.electronMaker.MVACuts = electronMVAcuts
process.photonMaker.MVACuts = photonMVAcuts

#for idmod in my_eleid_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
#for idmod in my_phoid_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


##########################
# Lepton post-processing #
##########################
# Muons
# Ghost cleaning
# See https://indico.cern.ch/event/203424/contributions/1489038/attachments/307148/428853/gpetrucc-ghosts-090812.pdf
process.cleanedMuons = cms.EDProducer(
   "PATMuonCleanerBySegments",
   src = cms.InputTag("slimmedMuons"),
   preselection = cms.string("track.isNonnull"), # Veto if this is not satisfied
   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"), # Avoid veto if these are satisfied
   fractionOfSharedSegments = cms.double(0.499)
   )
process.muonMaker.muonsInputTag = cms.InputTag("cleanedMuons")


#############
# Sequences ################################################################################################################################################
#############

# Gen. sequence
if not opts.data:
   process.genMakerSeq = cms.Sequence( process.genMaker )

# e/gamma sequence
process.egammaMakerSeq = cms.Sequence( process.heepIDVarValueMaps * process.egammaPostRecoSeq * process.electronMaker * process.photonMaker )

# Muon sequence
if opts.applyMuoncorr:
   process.correctedMuons.identifier = cms.string("RoccoR{}".format(opts.year))
   process.correctedMuons.isMC = cms.bool((not opts.data))

   process.cleanedMuons.src = cms.InputTag("correctedMuons")

   process.muonMakerSeq = cms.Sequence( process.correctedMuons * process.cleanedMuons * process.muonMaker )
else:
   process.muonMakerSeq = cms.Sequence( process.cleanedMuons * process.muonMaker )

# PF jets
if jecVersion != "":
   # Input to pfJetMaker is already fixed above
   process.pfJetMakerSeq = cms.Sequence( process.pileupJetIdUpdated + process.slimmedJetCorrFactors * process.slimmedCorrectedJets * process.pfJetMaker )
   process.pfJetPUPPIMakerSeq = cms.Sequence( process.slimmedPuppiJetCorrFactors * process.slimmedCorrectedJetsPuppi * process.pfJetPUPPIMaker )
   process.subJetMakerSeq = cms.Sequence( process.slimmedJetAK8CorrFactors * process.slimmedCorrectedJetsAK8 * process.subJetMaker )
else:
   process.pfJetMakerSeq = cms.Sequence( process.pileupJetIdUpdated + process.pfJetMaker )
   process.pfJetPUPPIMakerSeq = cms.Sequence( process.pfJetPUPPIMaker )
   process.subJetMakerSeq = cms.Sequence( process.subJetMaker )

# MET filter
if not opts.is80x:
   process.metFilterMakerSeq = cms.Sequence( process.ecalBadCalibReducedMINIAODFilter * process.metFilterMaker )
else:
   process.metFilterMakerSeq = cms.Sequence( process.metFilterMaker )

# PF MET
if hasattr(process, "fullPatMetSequenceModifiedMET"):
   process.pfmetMakerSeq = cms.Sequence( process.fullPatMetSequenceModifiedMET * process.pfmetMaker )
else:
   process.pfmetMakerSeq = cms.Sequence( process.pfmetMaker )
process.PFMETPath = cms.Path(process.pfmetMakerSeq) # Make as separate path

# Puppi MET
if hasattr(process, "fullPatMetSequenceModifiedPuppiMET"):
   process.pfmetpuppiMakerSeq = cms.Sequence( process.puppiMETSequence * process.fullPatMetSequenceModifiedPuppiMET * process.pfmetpuppiMaker )
else:
   process.pfmetpuppiMakerSeq = cms.Sequence( process.pfmetpuppiMaker )
process.PFMETPuppiPath = cms.Path(process.pfmetpuppiMakerSeq) # Make as separate path

# Isotracks
if opts.is80x:
   process.isoTrackMakerSeq = cms.Sequence( process.isoTrackMaker )
else:
   process.isoTrackMakerSeq = None



###################
# Build the paths ##########################################################################################################################################
###################

# steal some logic from https://github.com/cms-sw/cmssw/blob/CMSSW_10_4_X/PhysicsTools/NanoAOD/python/nano_cff.py
producers = [
        #process.eventMaker,
        process.lumiFilter, # filter after eventmaker so we get run lumi event branches at least, and event counts match up
        #process.muToTrigAssMaker if opts.triginfo else None,
        #process.elToTrigAssMaker if opts.triginfo else None,
        process.hltMaker if not opts.fastsim else None,
        #process.miniAODrhoSequence,
        process.metFilterMaker,
        process.vertexMaker,
        #process.secondaryVertexMaker,
        # process.pfCandidateMaker,
        process.muonMaker,
        process.electronMaker,
        process.pfJetMaker,
        #process.pfJetPUPPIMaker,
        process.subJetMaker,
        #process.pfmetMaker,
        #process.pfmetpuppiMaker,
        #process.photonMaker,
        #process.pftauMaker,
        process.isoTrackMaker,
        # Disable these temporarily until they are renewed
        process.genMaker if not opts.data else None,
        process.genJetMaker if not opts.data else None,
        #process.candToGenAssMaker if not opts.data else None,
        #process.pdfinfoMaker if not opts.data else None,
        #process.puSummaryInfoMaker if not opts.data else None,
        #process.hypDilepMaker,
        #process.sParmMaker if (opts.fastsim or opts.sparminfo) else None,
        ]

if opts.genxsecanalyzer and not opts.data:
    process.genxsecanalyzer = cms.EDAnalyzer("GenXSecAnalyzer")
    producers = [process.genxsecanalyzer]

total_path = None
for ip,producer in enumerate(producers):
   if producer is None: continue

   # Before adding the producers to the path, set the common attributes
   if hasattr(producer, "year"):
      setattr(producer, "year", cms.int32(opts.year))

   if ip == 0:
      total_path = producer
      continue

   if producer == process.isoTrackMaker:
      if not opts.is80x:
         total_path *= process.isoTrackMakerSeq
      continue

   if producer == process.metFilterMaker:
      if not (opts.fastsim and opts.is80x):
         total_path *= process.metFilterMakerSeq
      continue

   if producer == process.pfmetMaker:
      total_path *= process.pfmetMakerSeq
      continue

   if producer == process.pfmetpuppiMaker:
      total_path *= process.pfmetpuppiMakerSeq
      continue

   if producer == process.genMaker:
      if not opts.data:
         total_path *= process.genMakerSeq
      continue

   if producer == process.electronMaker:
      total_path *= process.egammaMakerSeq
      continue

   if producer == process.muonMaker:
      total_path *= process.muonMakerSeq
      continue

   if producer == process.pfJetMaker:
      total_path *= process.pfJetMakerSeq
      continue

   if producer == process.pfJetPUPPIMaker:
      total_path *= process.pfJetPUPPIMakerSeq
      continue

   if producer == process.subJetMaker:
      total_path *= process.subJetMakerSeq
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
# process.ProfilerService = cms.Service(
#    "ProfilerService",
#    firstEvent = cms.untracked.int32(2),
#    lastEvent = cms.untracked.int32(10),
#    paths = cms.untracked.vstring('p1')
#    )

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
#process.eventMaker.CMS3tag = cms.string('SUPPLY_CMS3_TAG')
#process.eventMaker.datasetName = cms.string('SUPPLY_DATASETNAME')

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
   process.cms3ntuple.is80X = cms.bool(opts.is80x)
   process.cms3ntuple.prefiringWeightsTag = cms.untracked.string(prefiringWeightsTag)
   process.cms3ntuple.keepGenParticles = cms.untracked.string(opts.keepGenParticles)
   process.cms3ntuple.keepGenJets = cms.bool(opts.keepGenJets)
   process.outpath = cms.EndPath(process.cms3ntuple)


if opts.dumpProcess:
   fprocdump = open(opts.output.replace('.root','_run_cfg.py'),'w')
   fprocdump.write(process.dumpPython())
   fprocdump.close()
