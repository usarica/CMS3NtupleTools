import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5JetID_cfi import ak5JetID

prunedUncorrectedCMS2Jets = cms.EDFilter("JetCollectionPruner",
   inputUncorrectedJetCollection = cms.InputTag("ak5CaloJets"),
   CaloJetCorrectorL2L3          = cms.string('ak5CaloL2L3'),
   uncorrectedJetPtCut           = cms.double(5.0), ##cut on the UNCORRECTED reco jets!!!!!
   usecorrectedCut               = cms.bool(False), #if you want to keep OR of passing uncorrected and corrected cuts
   correctedJetPtCut             = cms.double(10.0), ##cut on the CORRECTED reco jets!!!!!										 
)


cms2ak5JetID = ak5JetID.clone()
cms2ak5JetID.src = cms.InputTag("prunedUncorrectedCMS2Jets")
