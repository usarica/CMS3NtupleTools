import FWCore.ParameterSet.Config as cms

pftauMaker = cms.EDProducer("PFTauMaker",
	aliasPrefix = cms.untracked.string("taus_pf"),
        minleadPFChargedHadrCandPt = cms.double(5.),
        # PFTau collection
        pftausInputTag = cms.InputTag("fixedConePFTauProducer"),
      
)


