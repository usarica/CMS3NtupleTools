
import FWCore.ParameterSet.Config as cms
from TrackingTools.TrackAssociator.default_cfi import *

# see e.g.
# RecoEgamma/EgammaIsolationAlgos/python/electronTrackIsolationScone_cfi.py

elTkJuraIsoMaker = cms.EDProducer(
	"ElTkJuraIsoMaker",
	aliasPrefix = cms.untracked.string("els"),

	# input collections
	elsInputTag = cms.InputTag("electronMaker"),
	trackInputTag = cms.InputTag("trackMaker"),

	# isolation parameters
	trackIsoExtRadius   	= cms.double(0.3),

	# to validate with els_tkIso, set these 
	# to the default values
       trackIsoInRadius        = cms.double(0.015),
	trackIsoJurassicWidth 	= cms.double(0.01),
#	trackIsoJurassicWidth   = cms.double(0.0),
#        trackIsoInRadius        = cms.double(0.04),

	trackIsoMinPt       	= cms.double(0.7),
	trackIsoMind0       	= cms.double(9999.0),
	trackIsoMinz0       	= cms.double(0.2)
)



