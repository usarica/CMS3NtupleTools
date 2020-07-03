import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection,setupBTagging
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_81x, _chsalgos_94x, _chsalgos_102x

def configureJetCorrections(
   process,
   jetsrc,
   JECpayload,
   year,
   is80x,
   isData,
   updateJEC=True,
   updatePUJetID=False,
   updateBtagging=False
   ):
      if not JECpayload:
         raise RuntimeError("configureJetCorrections: JEC payload must be defined")

      currentJetSrc = jetsrc

      jetFinalSeq = "selectedFinalJetsSeq"+JECpayload
      setattr(process, jetFinalSeq, cms.Sequence())
      _process_jetFinalSeq = getattr(process, jetFinalSeq)

      pujetidnames = []

      ## Pile-up jet id
      if updatePUJetID:
         process.load("RecoJets.JetProducers.PileupJetID_cfi")
         # Updated PU jet id
         if is80x or year==2017 or year==2018:
            pujetidmod = "pileupJetIdUpdated"+JECpayload
            setattr(process, pujetidmod, process.pileupJetId.clone(
                  jets = cms.InputTag(currentJetSrc),
                  vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  inputIsCorrected = False,
                  applyJec = True,
                  jec = cms.string(JECpayload)
                  )
               )
            _process_pujetidmod = getattr(process, pujetidmod)
            if year == 2016:
               _process_pujetidmod.algos = cms.VPSet(_chsalgos_81x)
            elif year == 2017:
               _process_pujetidmod.algos = cms.VPSet(_chsalgos_94x)
            elif year == 2018:
               _process_pujetidmod.algos = cms.VPSet(_chsalgos_102x)
            if hasattr(_process_pujetidmod, "usePuppi"):
               _process_pujetidmod.usePuppi = cms.bool(("puppi" in JECpayload.lower()))
            elif "puppi" in JECpayload.lower():
               raise RuntimeError("Cannot apply PU jet id on Puppi jets in this CMSSW release")
            _process_jetFinalSeq *= _process_pujetidmod
            pujetidnames.append(pujetidmod)
         # Default PU jet id (recalculated for the JEC updates)
         pujetidmod = "pileupJetIdUpdated"+JECpayload+"Default"
         setattr(process, pujetidmod, process.pileupJetId.clone(
               jets = cms.InputTag(currentJetSrc),
               vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
               inputIsCorrected = False,
               applyJec = True,
               jec = cms.string(JECpayload)
               )
            )
         _process_pujetidmod_default = getattr(process, pujetidmod)
         if year >= 2016 and year<=2018:
            _process_pujetidmod_default.algos = cms.VPSet(_chsalgos_81x)
         if hasattr(_process_pujetidmod_default, "usePuppi"):
            _process_pujetidmod_default.usePuppi = cms.bool(("puppi" in JECpayload.lower()))
         elif "puppi" in JECpayload.lower():
            raise RuntimeError("Cannot apply PU jet id on Puppi jets in this CMSSW release")
         _process_jetFinalSeq *= _process_pujetidmod_default
         pujetidnames.append(pujetidmod)

      if updateBtagging:
         # First build a jet source that reverts JECs
         jetsrc_undojec_corrfactorproducer = currentJetSrc+"UndoJECCorrFactors"
         setattr(process, jetsrc_undojec_corrfactorproducer, updatedPatJetCorrFactors.clone(
               src = cms.InputTag(jetsrc),
               primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
               levels = cms.vstring(),
               payload = cms.string(JECpayload)
               )
            )
         _process_jetsrc_undojec_corrfactorproducer = getattr(process, jetsrc_undojec_corrfactorproducer)

         jetsrc_undojec = currentJetSrc+"UndoJEC"
         setattr(process, jetsrc_undojec, updatedPatJets.clone(
               jetSource = cms.InputTag(currentJetSrc),
               addJetCorrFactors = cms.bool(True),
               jetCorrFactorsSource = cms.VInputTag(cms.InputTag(jetsrc_undojec_corrfactorproducer)),
               addBTagInfo = cms.bool(True),
               addDiscriminators = cms.bool(True),
               addTagInfos = cms.bool(False),
               discriminatorSources = cms.VInputTag(),
               tagInfoSources = cms.VInputTag()
               )
            )
         _process_jetsrc_undojec = getattr(process, jetsrc_undojec)

         if updatePUJetID:
            for pujetidname in pujetidnames:
               _process_jetsrc_undojec.userData.userFloats.src += [pujetidname+":fullDiscriminant"]
               _process_jetsrc_undojec.userData.userInts.src += [pujetidname+":fullId"]

         jetsrc_finaljec_corrfactorproducer = currentJetSrc+"FinalJECCorrFactors"
         setattr(process, jetsrc_finaljec_corrfactorproducer, updatedPatJetCorrFactors.clone(
               src = cms.InputTag(jetsrc_undojec),
               primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
               levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute']),
               payload = cms.string(JECpayload)
               )
            )
         _process_jetsrc_finaljec_corrfactorproducer = getattr(process, jetsrc_finaljec_corrfactorproducer)
         ## Data applies L2L3Residual corrections as well
         if isData:
            _process_jetsrc_finaljec_corrfactorproducer.levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])

         jetsrc_finaljets = currentJetSrc+"FinalJets"
         setattr(process, jetsrc_finaljets, updatedPatJets.clone(
               jetSource = cms.InputTag(jetsrc_undojec),
               addJetCorrFactors = cms.bool(True),
               jetCorrFactorsSource = cms.VInputTag(cms.InputTag(jetsrc_finaljec_corrfactorproducer)),
               addBTagInfo = cms.bool(True),
               addDiscriminators = cms.bool(True),
               addTagInfos = cms.bool(False),
               discriminatorSources = cms.VInputTag(),
               tagInfoSources = cms.VInputTag()
               )
            )
         _process_jetsrc_finaljets = getattr(process, jetsrc_finaljets)

         deep_discriminators = []
         if is80x:
            deep_discriminators.extend(['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfDeepCSVJetTags:probc'])
         deep_discriminators.extend([
            'pfDeepFlavourJetTags:probb',
            'pfDeepFlavourJetTags:probbb',
            'pfDeepFlavourJetTags:problepb',
            'pfDeepFlavourJetTags:probc',
            'pfDeepFlavourJetTags:probuds',
            'pfDeepFlavourJetTags:probg'
            ])
         setupBTagging(
            process,
            jetSource = cms.InputTag(jetsrc_undojec),
            pfCandidates = cms.InputTag('packedPFCandidates'),
            explicitJTA = False,
            pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
            svSource = cms.InputTag('slimmedSecondaryVertices'),
            elSource = cms.InputTag('slimmedElectrons'),
            muSource = cms.InputTag('slimmedMuons'),
            runIVF = False,
            tightBTagNTkHits = False,
            loadStdRecoBTag = False,
            svClustering = False,
            fatJets = cms.InputTag(''),
            groomedFatJets = cms.InputTag(''),
            algo = 'AK',
            rParam = 0.4,
            btagDiscriminators = deep_discriminators,
            btagInfos = [],
            patJets = _process_jetsrc_finaljets,
            labelName = '',
            btagPrefix = '',
            postfix = JECpayload
            )

         _process_jetFinalSeq *= _process_jetsrc_undojec_corrfactorproducer
         _process_jetFinalSeq *= _process_jetsrc_undojec
         _process_jetFinalSeq *= getattr(process, "pfImpactParameterTagInfos"+JECpayload)
         _process_jetFinalSeq *= getattr(process, "pfInclusiveSecondaryVertexFinderTagInfos"+JECpayload)
         _process_jetFinalSeq *= getattr(process, "pfDeepCSVTagInfos"+JECpayload)
         if is80x:
            _process_jetFinalSeq *= getattr(process, "pfDeepCSVJetTags"+JECpayload)
         _process_jetFinalSeq *= getattr(process, "pfDeepFlavourTagInfos"+JECpayload)
         _process_jetFinalSeq *= getattr(process, "pfDeepFlavourJetTags"+JECpayload)
         _process_jetFinalSeq *= _process_jetsrc_finaljec_corrfactorproducer
         _process_jetFinalSeq *= _process_jetsrc_finaljets
         currentJetSrc = jetsrc_finaljets
      else:
         jetsrc_finaljec_corrfactorproducer = currentJetSrc+"FinalJECCorrFactors"
         setattr(process, jetsrc_finaljec_corrfactorproducer, updatedPatJetCorrFactors.clone(
               src = cms.InputTag(currentJetSrc),
               primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
               levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute']),
               payload = cms.string(JECpayload)
               )
            )
         _process_jetsrc_finaljec_corrfactorproducer = getattr(process, jetsrc_finaljec_corrfactorproducer)
         ## Data applies L2L3Residual corrections as well
         if isData:
            _process_jetsrc_finaljec_corrfactorproducer.levels = cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])

         jetsrc_finaljets = currentJetSrc+"FinalJets"
         setattr(process, jetsrc_finaljets, updatedPatJets.clone(
               jetSource = cms.InputTag(currentJetSrc),
               addJetCorrFactors = cms.bool(True),
               jetCorrFactorsSource = cms.VInputTag(cms.InputTag(jetsrc_finaljec_corrfactorproducer)),
               addBTagInfo = cms.bool(True),
               addDiscriminators = cms.bool(True),
               addTagInfos = cms.bool(False),
               discriminatorSources = cms.VInputTag(),
               tagInfoSources = cms.VInputTag()
               )
            )
         _process_jetsrc_finaljets = getattr(process, jetsrc_finaljets)

         if updatePUJetID:
            for pujetidname in pujetidnames:
               _process_jetsrc_finaljets.userData.userFloats.src += [pujetidname+":fullDiscriminant"]
               _process_jetsrc_finaljets.userData.userInts.src += [pujetidname+":fullId"]

         _process_jetFinalSeq *= _process_jetsrc_finaljec_corrfactorproducer
         _process_jetFinalSeq *= _process_jetsrc_finaljets
         currentJetSrc = jetsrc_finaljets


      jetsrc_selectedfinaljets = "selectedFinalJets"+JECpayload
      setattr(process, jetsrc_selectedfinaljets, cms.EDFilter(
            "PATJetSelector",
            cut = cms.string(''),
            cutLoose = cms.string(''),
            nLoose = cms.uint32(0),
            src = cms.InputTag(currentJetSrc)
            )
         )
      _process_jetsrc_selectedfinaljets = getattr(process, jetsrc_selectedfinaljets)
      _process_jetFinalSeq *= _process_jetsrc_selectedfinaljets

      print "configureJetCorrections: Final jet sequence for {} jets:".format(JECpayload)
      print _process_jetFinalSeq

      return _process_jetFinalSeq
