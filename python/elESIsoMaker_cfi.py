import FWCore.ParameterSet.Config as cms

elESIsoMaker = cms.EDProducer("ElESIsoMaker",
	aliasPrefix = cms.untracked.string("els"),

  	electronsInputTag 	= cms.InputTag("gsfElectrons"),
        esHitsInputTag 		= cms.InputTag("ecalPreshowerRecHit:EcalRecHitsES"),
	intRadius = cms.double(3.0), # crystals!
	etaSlice = cms.double(1.5),  # crystals! NOTE- surely that is too small!!!
	useNumCrystals = cms.bool(True)

)

