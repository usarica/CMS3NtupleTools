import FWCore.ParameterSet.Config as cms

calotauMaker = cms.EDFilter("CaloTauMaker",
	aliasPrefix = cms.untracked.string("taus_calo"),
        minleadTrackPt = cms.double(5.),
        # CaloTau collection
        calotausInputTag = cms.InputTag("caloRecoTauProducer")
)


