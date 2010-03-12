import FWCore.ParameterSet.Config as cms

pftauMaker = cms.EDFilter("PFTauMaker",
	aliasPrefix = cms.untracked.string("taus_pf"),
        # PFTau collection
        pftausInputTag = cms.InputTag("fixedConePFTauProducer"),
      
)


