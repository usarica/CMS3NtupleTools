import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5JetID_cfi import ak5JetID

prunedUncorrectedCMS2Jets = cms.EDProducer("JetCollectionPruner",
                                         inputUncorrectedJPTJetCollection = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
                                         inputUncorrectedPFJetCollection  = cms.InputTag("ak5PFJets"),
                                         inputUncorrectedTrkJetCollection = cms.InputTag("ak5TrackJets"),
                                         uncorrectedJPTJetPtCut           = cms.double(8.0), ##cut on uncorrected JPT jets!!!!!
                                         uncorrectedTrkJetPtCut           = cms.double(5.0),                                                                                   
                                         uncorrectedPFJetPtCut            = cms.double(7.0)
)
                                         

cms2ak5JetID = ak5JetID.clone()
cms2ak5JetID.src = cms.InputTag("prunedUncorrectedCMS2Jets", "calojet")
