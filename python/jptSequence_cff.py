import FWCore.ParameterSet.Config as cms

from RecoJets.JetPlusTracks.JetPlusTrackCorrections_cff import *

from CMS2.NtupleMaker.jptMaker_cfi import jptMaker

JetPlusTrackZSPCorJetAntiKt5.src = cms.InputTag("prunedUncorrectedCMS2Jets")

JPTCorrections = cms.Sequence(JetPlusTrackCorrectionsAntiKt5 * jptMaker)
