import FWCore.ParameterSet.Config as cms

from RecoJets.JetAssociationProducers.trackExtrapolator_cfi import *
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("prunedUncorrectedCMS2Jets")
ak5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("prunedUncorrectedCMS2Jets")
ak5JetExtender.jets = cms.InputTag("prunedUncorrectedCMS2Jets")

from RecoJets.JetPlusTracks.JetPlusTrackCorrections_cff import *
JetPlusTrackZSPCorJetAntiKt5.src = cms.InputTag("prunedUncorrectedCMS2Jets")

from CMS2.NtupleMaker.jptMaker_cfi import jptMaker

JPTCorrections = cms.Sequence(trackExtrapolator * ak5JTA * JetPlusTrackCorrectionsAntiKt5 * jptMaker)
