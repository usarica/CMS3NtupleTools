import FWCore.ParameterSet.Config as cms

calotauMaker = cms.EDFilter("CaloTauMaker",
	aliasPrefix = cms.untracked.string("taus_calo"),
      
        # CaloTau collection
        calotausInputTag = cms.InputTag("caloRecoTauProducer")
)


