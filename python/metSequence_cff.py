import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff import *
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorIC5CaloJet

metMuonJESCorSC5CMS2                     = metJESCorIC5CaloJet.clone()
metMuonJESCorSC5CMS2.inputUncorJetsLabel = "prunedUncorrectedCMS2Jets"
metMuonJESCorSC5CMS2.corrector           = "L2L3JetCorrectorSC5Calo"
metMuonJESCorSC5CMS2.inputUncorMetLabel  = "corMetGlobalMuons"

metCorSequence = cms.Sequence(metMuonJESCorSC5CMS2)

