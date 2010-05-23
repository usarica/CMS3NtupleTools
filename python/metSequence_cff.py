import FWCore.ParameterSet.Config as cms


from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorIC5CaloJet
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

metMuonJESCorAK5CMS2                     = metJESCorIC5CaloJet.clone()
metMuonJESCorAK5CMS2.inputUncorJetsLabel = "prunedUncorrectedCMS2Jets"
metMuonJESCorAK5CMS2.corrector           = "ak5CaloL2L3"
metMuonJESCorAK5CMS2.inputUncorMetLabel  = "corMetGlobalMuons"

from RecoJets.JetProducers. CaloTowerSchemeB_cfi import *
from RecoMET.METProducers.CaloMET_cfi import met

cms2towerMaker = towerMaker.clone()
metNew = met.clone()
met.src = cms.InputTag("cms2towerMaker")

from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *
cms2HBHENoiseFilterResultProducer = HBHENoiseFilterResultProducer.clone()

metCorSequence = cms.Sequence(metMuonJESCorAK5CMS2 * cms2towerMaker * metNew * cms2HBHENoiseFilterResultProducer)
