import FWCore.ParameterSet.Config as cms

calotauMaker = cms.EDFilter("CaloTauMaker",
      
        # CaloTau collection
        calotausInputTag = cms.InputTag("caloRecoTauProducer")
)


