import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.jetCollectionPruner_cfi import prunedUncorrectedCMS2Jets

prunedUncorrectedCMS2scJets = prunedUncorrectedCMS2Jets.clone()
prunedUncorrectedCMS2scJets.inputUncorrectedJetCollection = cms.InputTag("sisCone5CaloJets")

cms2scCaloJetSequence = cms.Sequence(prunedUncorrectedCMS2scJets)
