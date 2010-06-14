# The following comments couldn't be translated into the new config version:

# more, if you've got these in the path:

import FWCore.ParameterSet.Config as cms

hypFilter = cms.EDProducer("TheFilter",
    # this syntax means: keep events that contain more than 0
    # entries in the collection hypDilepMaker:hyptype  
    filterExpressions = cms.VInputTag(cms.InputTag("hypDilepMaker","hyptype"), cms.InputTag("hypTrilepMaker","hyptrilepfirsttype"), cms.InputTag("hypQuadlepMaker","hypquadlepfirsttype"))
)


