import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.ZSPJetCorrections219_cff import *
from JetMETCorrections.Configuration.JetPlusTrackCorrections_cff import *
from CMS2.NtupleMaker.jptMaker_cfi import jptMaker

ZSPJetCorJetAntiKt5.src = cms.InputTag("prunedUncorrectedCMS2Jets")

JPTCorrections = cms.Sequence(ZSPJetCorrectionsAntiKt5 * JetPlusTrackCorrectionsAntiKt5 * jptMaker)
