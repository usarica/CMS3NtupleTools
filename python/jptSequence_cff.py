import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.ZSPJetCorrections219_cff import *
from JetMETCorrections.Configuration.JetPlusTrackCorrections_cff import *
from JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff import *
from CMS2.NtupleMaker.jptMaker_cfi import jptMaker

ZSPJetCorJetIcone5.src = cms.InputTag("prunedUncorrectedCMS2Jets")

JPTCorrections = cms.Sequence(ZSPJetCorrections * JetPlusTrackCorrections * L2L3CorJetIC5JPT * jptMaker)
