import FWCore.ParameterSet.Config as cms

pftauMaker = cms.EDFilter("PFTauMaker",
        # PFTau collection
        pftausInputTag = cms.InputTag("fixedConePFTauProducer"),
      
)


