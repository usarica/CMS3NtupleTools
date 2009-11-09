import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.L2L3Corrections_Summer09_cff import *
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorIC5CaloJet

metMuonJESCorAK5CMS2                     = metJESCorIC5CaloJet.clone()
metMuonJESCorAK5CMS2.inputUncorJetsLabel = "prunedUncorrectedCMS2Jets"
metMuonJESCorAK5CMS2.corrector           = "L2L3JetCorrectorAK5Calo"
metMuonJESCorAK5CMS2.inputUncorMetLabel  = "corMetGlobalMuons"

metCorSequence = cms.Sequence(metMuonJESCorAK5CMS2)

