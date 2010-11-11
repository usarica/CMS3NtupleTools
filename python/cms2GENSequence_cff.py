#only CMS2 gen makers

import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.candToGenAssMaker_cfi     import *
from CMS2.NtupleMaker.flavorHistorySequence_cfi import *
from CMS2.NtupleMaker.genJetSequence_cff        import *
from CMS2.NtupleMaker.genMaker_cfi              import *
from CMS2.NtupleMaker.hypGenMaker_cfi           import *
from CMS2.NtupleMaker.pdfinfoMaker_cfi          import *
from CMS2.NtupleMaker.genJetMaker_cfi           import *
from CMS2.NtupleMaker.puSummaryInfoMaker_cfi    import *

#This should be in the config
from CMS2.NtupleMaker.dilepGenFilter_cfi import dilepGenFilter


cms2GENSequence     = cms.Sequence(genMaker * pdfinfoMaker * 
                                   genJetSequence * CMS2FlavorHistorySequence * 
                                   candToGenAssMaker * genJetMaker * 
                                   hypGenMaker * puSummaryInfoMaker)
