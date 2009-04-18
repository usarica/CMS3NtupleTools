import FWCore.ParameterSet.Config as cms

# File: JetCorrections.cff
# Author: O. Kodolova
# Date: 1/24/07
#
# Jet corrections with tracks for the icone5 jets with ZSP corrections.
# ... adapted by DLE for SisCone jets ...
# 

from JetMETCorrections.Configuration.JetCorrectionsRecord_cfi import *
from RecoJets.Configuration.RecoJetAssociations_cff import *

JetPlusTrackZSPCorrectorSisCone = cms.ESSource("JetPlusTrackCorrectionService",
    JetTrackCollectionAtCalo = cms.InputTag("ZSPsisConeJetTracksAssociatorAtCaloFace"),
    respalgo = cms.int32(5),
    JetTrackCollectionAtVertex = cms.InputTag("ZSPsisConeJetTracksAssociatorAtVertex"),
    muonSrc = cms.InputTag("muons"),
    AddOutOfConeTracks = cms.bool(True),
    NonEfficiencyFile = cms.string('CMSSW_167_TrackNonEff'),
    NonEfficiencyFileResp = cms.string('CMSSW_167_TrackLeakage'),
    ResponseFile = cms.string('CMSSW_167_response'),
    label = cms.string('JetPlusTrackZSPCorrectorSisCone'),
    TrackQuality = cms.string('highPurity'),
    UseQuality = cms.bool(True)
)

JetPlusTrackZSPCorJetSisCone = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("ZSPJetCorJetSisCone"),
    correctors = cms.vstring('JetPlusTrackZSPCorrectorSisCone'),
    alias = cms.untracked.string('JetPlusTrackZSPCorJetSisCone')
)

from RecoJets.JetAssociationProducers.iterativeCone5JTA_cff import*

ZSPsisConeJetTracksAssociatorAtVertex = iterativeCone5JetTracksAssociatorAtVertex.clone()
ZSPsisConeJetTracksAssociatorAtVertex.jets = cms.InputTag("ZSPJetCorJetSisCone")

ZSPsisConeJetTracksAssociatorAtCaloFace = iterativeCone5JetTracksAssociatorAtCaloFace.clone()
ZSPsisConeJetTracksAssociatorAtCaloFace.jets = cms.InputTag("ZSPJetCorJetSisCone")

ZSPsisConeJetExtender = iterativeCone5JetExtender.clone()
ZSPsisConeJetExtender.jets = cms.InputTag("ZSPJetCorJetSisCone")
ZSPsisConeJetExtender.jet2TracksAtCALO = cms.InputTag("ZSPsisConeJetTracksAssociatorAtCaloFace")
ZSPsisConeJetExtender.jet2TracksAtVX = cms.InputTag("ZSPsisConeJetTracksAssociatorAtVertex")


ZSPrecoJetAssociationsSisCone = cms.Sequence(ZSPsisConeJetTracksAssociatorAtVertex*ZSPsisConeJetTracksAssociatorAtCaloFace*ZSPsisConeJetExtender)

JetPlusTrackCorrectionsSisCone = cms.Sequence(ZSPrecoJetAssociationsSisCone*JetPlusTrackZSPCorJetSisCone)

