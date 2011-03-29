# The following comments couldn't be translated into the new config version:

# more, if you've got these in the path:

import FWCore.ParameterSet.Config as cms

hypDiLeptonFilter = cms.EDFilter("HypDilepPruner",
    hyp_lt_p4       = cms.InputTag("hypDilepMaker","hypltp4"),
    hyp_ll_p4       = cms.InputTag("hypDilepMaker","hypllp4"),
    hyp_sumJetPt    = cms.InputTag("hypDilepMaker","hypsumJetPt"),
    minNominalLTPt  = cms.double(20.),
    minNominalLLPt  = cms.double(10.),
    minSusyLTPt     = cms.double(5.),
    minSusyLLPt     = cms.double(5.),
    minSusySumJetPt = cms.double(100.)
)

hypOtherFilter = cms.EDFilter("TheFilter",
    # this syntax means: keep events that contain more than 0
    # entries in the collection hypDilepMaker:hyptype  
    filterExpressions = cms.VInputTag(cms.InputTag("hypTrilepMaker","hyptrilepfirsttype"),
                                      cms.InputTag("hypQuadlepMaker","hypquadlepfirsttype"))
)
