# The following comments couldn't be translated into the new config version:

# -*-sh-*-

import FWCore.ParameterSet.Config as cms

process = cms.Process("CMS2")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cff")
process.l1extraParticles.muonSource = cms.InputTag("gtDigis")

process.GlobalTag.globaltag = "STARTUP_V4::All"

process.load("PhysicsTools.PatAlgos.patLayer0_cff")

process.load("PhysicsTools.PatAlgos.patLayer1_cff")



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
    input = cms.untracked.int32(100)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_2_1_9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/1A0FD639-1B86-DD11-A3C0-000423D99614.root')
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
    fileName = cms.untracked.string('ntuple_test.root')
)

process.JetCorrectionsExtra = cms.Sequence(process.L4EMFJetCorJetIcone5*process.MCJetCorJetIcone5)
process.JetCorrection = cms.Sequence(process.JetCorrectionsExtra)
process.pat = cms.Sequence(process.patLayer0*process.patLayer1)
process.makers = cms.Sequence(process.beamSpotMaker*process.muonMaker*process.electronMaker*process.jetMaker*process.trackMaker)
process.patmakers = cms.Sequence(process.patMuonMaker*process.patElectronMaker*process.patJetMaker)
process.assmakers = cms.Sequence(process.jetToMuAssMaker*process.jetToElAssMaker*process.muToElsAssMaker*process.candToGenAssMaker*process.muToJetAssMaker*process.muToTrackAssMaker*process.elToTrackAssMaker*process.elToMuAssMaker*process.trackToMuonAssMaker*process.trackToElsAssMaker)
process.generalmakers = cms.Sequence(process.l1extraParticles*process.eventMaker*process.metMaker*process.l1DigiMaker*process.genMaker)
process.hypmaker = cms.Sequence(process.hypTrilepMaker*process.hypDilepMaker*process.hypQuadlepMaker)
process.cms2 = cms.Sequence(process.generalmakers*process.makers*process.patmakers*process.assmakers*process.hypmaker)
process.p = cms.Path(process.JetCorrection*process.pat*process.cms2)
process.outpath = cms.EndPath(process.outMod)

