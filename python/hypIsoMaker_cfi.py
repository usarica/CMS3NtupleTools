import FWCore.ParameterSet.Config as cms
from TrackingTools.TrackAssociator.default_cfi import *

hypIsoMaker = cms.EDProducer(
	"HypIsoMaker",
	TrackAssociatorParameterBlock, #this is read as TrackAssociatorParameters (don't ask me why)

	#This is whether to use the egamma/muon isolations and subtract the isolation veto of the
	#other hypothesis (False), or recompute the isolation variable from hits (True)
	#recomputeEcalIso = cms.bool(False),
	#recomputeTrckIso = cms.bool(False),
	recomputeEcalIso = cms.bool(True),
	recomputeTrckIso = cms.bool(True),
	
	emObjectProducer = cms.InputTag("uniqueElectrons"),
    muonsInputTag = cms.InputTag("muons"),
	#cms2 collections
	cms2elsInputTag = cms.InputTag("electronMaker"),
	cms2musInputTag = cms.InputTag("muonMaker"),
	trackInputTag = cms.InputTag("trackMaker"),
	hypInputTag = cms.InputTag("hypDilepMaker"),
	caloTowersInputTag = cms.InputTag("towerMaker"),

	ecalBarrelRecHitProducer 	= cms.InputTag("ecalRecHit"),
	ecalBarrelRecHitCollection 	= cms.InputTag("EcalRecHitsEB"),
	ecalEndcapRecHitProducer 	= cms.InputTag("ecalRecHit"),
	ecalEndcapRecHitCollection 	= cms.InputTag("EcalRecHitsEE"),
	
	elsEcalVetoRadBarrel 	= cms.double(0.045), #els inner cone size--Barrel
	elsEcalVetoRadEndcap 	= cms.double(0.07),  #els inner cone size--Endcap
	elsEcalExtCone 		    = cms.double(0.4),   #els outer cone size		       
	musEcalVetoRadBarrel 	= cms.double(0.07),  #mus inner cone size--Barrel: note, larger b'c uses towers
	musEcalVetoRadEndcap 	= cms.double(0.07),  #mus inner cone size--Endcap: towers are larger in endcap(?)
	musEcalExtCone          = cms.double(0.3),   #mus outer cone size--use smaller veto for mus

	#electron thresholds--applied to rechits
	elsetMinBarrel 		= cms.double(-9999), #cut on energy
	elseMinBarrel 		= cms.double(0.08),  #cuts on et
	elsetMinEndcap 		= cms.double(-9999),
	elseMinEndcap 		= cms.double(0.3),  
	#elsetMinBarrel 		= cms.double(0.08), #cuts on et
	#elseMinBarrel 		= cms.double(-9999),
	#elsetMinEndcap 		= cms.double(0.3),
	#elseMinEndcap 		= cms.double(-9999),

	#muon thresholds--applied to towers
	#musetMinBarrel 		= cms.double(-9999), #cuts on energy
	#museMinBarrel 		= cms.double(0.2),
	#musetMinEndcap 		= cms.double(-9999),
	#museMinEndcap 		= cms.double(0.2),
	musetMinBarrel 		= cms.double(0.2),  #cuts on et
	museMinBarrel 		= cms.double(-9999),
	musetMinEndcap 		= cms.double(0.2),
	museMinEndcap 		= cms.double(-9999),
	
	jurassicWidth = cms.double(0.02), #dEta half-strip width
	
	useIsolEt = cms.bool(True),         #we always want et so we can compare to pt

	trackIsoExtRadius   = cms.double(0.3),
	trackIsoElsInRadius = cms.double(0.015),
	trackIsoMusInRadius = cms.double(0.01),
	trackIsoMinPt       = cms.double(1.0),
	trackIsoMind0       = cms.double(0.1), #mu defaults
	trackIsoMinz0       = cms.double(0.2),
	#el defaults are same--see RecoEgamma/EgammaIsolationAlgos/python/eleTrackExtractorBlocks_cff.py
	#trackIsoMind0       = cms.double(999.0), #my results are worse for no d0 than 0.1 d0... don't know why
	#trackIsoMinz0       = cms.double(999.0),
)
