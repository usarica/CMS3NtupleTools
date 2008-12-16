# The following comments couldn't be translated into the new config version:

# -*-sh-*-

import FWCore.ParameterSet.Config as cms

process = cms.Process("CMS2")

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.8 $'),
        annotation = cms.untracked.string('CMS2'),
        name = cms.untracked.string('CMS2 test configuration')
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cff")
#process.l1extraParticles.muonSource = cms.InputTag("gtDigis")

process.GlobalTag.globaltag = "IDEAL_V9::All"

#replace with the patched PAT layer 0 and layer1,
#patched for the electron isolation
process.load("CMS2.NtupleMaker.PATElecIsoPatch_cfi")

# CMS2 includes
process.load("CMS2.NtupleMaker.beamSpotMaker_cfi")

process.load("CMS2.NtupleMaker.eventMaker_cfi")

process.load("CMS2.NtupleMaker.jetMaker_cfi")

process.load("CMS2.NtupleMaker.jetToMuAssMaker_cfi")

process.load("CMS2.NtupleMaker.jetToElAssMaker_cfi")

process.load("CMS2.NtupleMaker.patJetMaker_cfi")

process.load("CMS2.NtupleMaker.electronMaker_cfi")

process.load("CMS2.NtupleMaker.patElectronMaker_cfi")

process.load("CMS2.NtupleMaker.hypTrilepMaker_cfi")

process.load("CMS2.NtupleMaker.metMaker_cfi")

process.load("CMS2.NtupleMaker.genMaker_cfi")

process.load("CMS2.NtupleMaker.candToGenAssMaker_cfi")

process.load("CMS2.NtupleMaker.muonMaker_cfi")

process.load("CMS2.NtupleMaker.trackMaker_cfi")

process.load("CMS2.NtupleMaker.patMETMaker_cfi")

process.load("CMS2.NtupleMaker.patMuonMaker_cfi")

process.load("CMS2.NtupleMaker.muToElsAssMaker_cfi")

process.load("CMS2.NtupleMaker.muToJetAssMaker_cfi")

process.load("CMS2.NtupleMaker.muToTrackAssMaker_cfi")

process.load("CMS2.NtupleMaker.elToTrackAssMaker_cfi")

process.load("CMS2.NtupleMaker.elToMuAssMaker_cfi")

process.load("CMS2.NtupleMaker.trackToMuonAssMaker_cfi")

process.load("CMS2.NtupleMaker.trackToElsAssMaker_cfi")

process.load("CMS2.NtupleMaker.hypDilepMaker_cfi")

process.load("CMS2.NtupleMaker.hypQuadlepMaker_cfi")

process.load("CMS2.NtupleMaker.theFilter_cfi")

#process.Timing = cms.Service("Timing")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring('file:///uscms_data/d1/slava77/tauola-16AAC418-218A-DD11-AC33-001F2908F0E4.root')
   #secondaryFileNames = cms.untracked.vstring('/store/mc/Summer08/WJets-madgraph/GEN-SIM-RAW/IDEAL_V9_v1/0030/D4CD0886-BA9E-DD11-8B40-003048770C6C.root',
   #                                           '/store/mc/Summer08/WJets-madgraph/GEN-SIM-RAW/IDEAL_V9_v1/0030/7EC5CAC6-1A9F-DD11-9811-003048770BAA.root',
   #                                           '/store/mc/Summer08/WJets-madgraph/GEN-SIM-RAW/IDEAL_V9_v1/0032/8A31B75C-F29F-DD11-A96E-0002B3E92671.root',
   #                                           '/store/mc/Summer08/WJets-madgraph/GEN-SIM-RAW/IDEAL_V9_v1/0030/0C79D330-169F-DD11-8585-003048770DBE.root')
)

process.MCJetCorrectorIcone5 = cms.ESSource("MCJetCorrectionService",
    tagName = cms.string('CMSSW_152_iterativeCone5'),
    label = cms.string('MCJetCorrectorIcone5')
)

process.L4EMFJetCorrectorIcone5 = cms.ESSource("L4EMFCorrectionService",
    tagName = cms.string('Spring07_L4EMF_Iterative_Cone_05'),
    label = cms.string('L4EMFJetCorrectorIcone5')
)

process.L4EMFJetCorrectorIcone7 = cms.ESSource("L4EMFCorrectionService",
    tagName = cms.string('Spring07_L4EMF_Iterative_Cone_05'),
    label = cms.string('L4EMFJetCorrectorIcone7')
)

process.L4EMFJetCorrectorMcone5 = cms.ESSource("L4EMFCorrectionService",
    tagName = cms.string('Spring07_L4EMF_Iterative_Cone_05'),
    label = cms.string('L4EMFJetCorrectorMcone5')
)

process.L4EMFJetCorrectorMcone7 = cms.ESSource("L4EMFCorrectionService",
    tagName = cms.string('Spring07_L4EMF_Iterative_Cone_05'),
    label = cms.string('L4EMFJetCorrectorMcone7')
)

process.L4EMFJetCorJetIcone5 = cms.EDProducer("JetCorrectionProducer",
    src = cms.InputTag("iterativeCone5CaloJets"),
    correctors = cms.vstring('MCJetCorrectorIcone5', 
        'L4EMFJetCorrectorIcone5'),
    alias = cms.untracked.string('L4EMFJetCorJetIcone5')
)

process.MCJetCorJetIcone5 = cms.EDProducer("JetCorrectionProducer",
    src = cms.InputTag("iterativeCone5CaloJets"),
    correctors = cms.vstring('MCJetCorrectorIcone5'),
    alias = cms.untracked.string('MCJetCorJetIcone5')
)


process.l1DigiMaker = cms.EDFilter("L1DigiMaker")

process.outMod = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *',
        'keep *_*Maker_*_CMS2'), 
    fileName = cms.untracked.string('ntuple_ttbar_nofilter.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('USER')
    )
)

process.triggerEventMaker = cms.EDProducer("TriggerEventMaker")

process.JetCorrectionsExtra = cms.Sequence(process.L4EMFJetCorJetIcone5*process.MCJetCorJetIcone5)
process.JetCorrection = cms.Sequence(process.JetCorrectionsExtra)
process.pat = cms.Sequence(process.patchPATSequence)
process.makers = cms.Sequence(process.beamSpotMaker*process.muonMaker*process.electronMaker*process.jetMaker*process.trackMaker)
process.patmakers = cms.Sequence(process.patMuonMaker*process.patElectronMaker*process.patJetMaker*process.patMETMaker)
process.assmakers = cms.Sequence(process.jetToMuAssMaker*process.jetToElAssMaker*process.muToElsAssMaker*process.candToGenAssMaker*process.muToJetAssMaker*process.muToTrackAssMaker*process.elToTrackAssMaker*process.elToMuAssMaker*process.trackToMuonAssMaker*process.trackToElsAssMaker)
process.trigprimmakers = cms.Sequence(process.l1DigiMaker*process.triggerEventMaker)
process.generalmakers = cms.Sequence(process.eventMaker*process.metMaker*process.genMaker)
process.hypmaker = cms.Sequence(process.hypTrilepMaker*process.hypDilepMaker*process.hypQuadlepMaker)
process.cms2 = cms.Sequence(process.generalmakers*process.trigprimmakers*process.makers*process.patmakers*process.assmakers*process.hypmaker)
process.p = cms.Path(process.JetCorrection*process.pat*process.cms2)
#process.outpath = cms.EndPath(process.theFilter*process.outMod)
process.outpath = cms.EndPath(process.outMod)

