import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.jetCollectionPruner_cfi import prunedUncorrectedCMS2Jets
from RecoJets.JetProducers.sc5JetID_cfi import sc5JetID


prunedUncorrectedCMS2scJets = prunedUncorrectedCMS2Jets.clone()
prunedUncorrectedCMS2scJets.inputUncorrectedJetCollection = cms.InputTag("sisCone5CaloJets")


cms2sc5JetID = sc5JetID.clone()
cms2sc5JetID.src = cms.InputTag("prunedUncorrectedCMS2scJets")


cms2scCaloJetSequence = cms.Sequence(prunedUncorrectedCMS2scJets*cms2sc5JetID)
