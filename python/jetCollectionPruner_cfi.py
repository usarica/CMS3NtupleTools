import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5JetID_cfi import ak5JetID

prunedUncorrectedCMS2Jets = cms.EDFilter("JetCollectionPruner",
   inputUncorrectedJetCollection = cms.InputTag("ak5CaloJets"),
   uncorrectedJetPtCut           = cms.double(5.0) ##cut on the UNCORRECTED reco jets!!!!!
)


cms2ak5JetID = ak5JetID.clone()
cms2ak5JetID.src = cms.InputTag("prunedUncorrectedCMS2Jets")
