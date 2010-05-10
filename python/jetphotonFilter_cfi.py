import FWCore.ParameterSet.Config as cms

jetphotonFilter = cms.EDFilter("TheFilter",
    # this syntax means: keep events that contain more than 0 entries in the collections listed below
    # branches must be ints or floats--no lorentzvectors allowed bc not implemented in TheFilter							   
    filterExpressions = cms.VInputTag(cms.InputTag("jetMaker","jetsemFrac"), cms.InputTag("photonMaker","photonscindex"))
)


