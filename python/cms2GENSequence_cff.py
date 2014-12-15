#only CMS2 gen makers

import FWCore.ParameterSet.Config as cms

from CMS3.NtupleMaker.candToGenAssMaker_cfi     import *
from CMS3.NtupleMaker.flavorHistorySequence_cfi import *
from CMS3.NtupleMaker.genJetSequence_cff        import *
from CMS3.NtupleMaker.genMaker_cfi              import *
from CMS3.NtupleMaker.hypGenMaker_cfi           import *
from CMS3.NtupleMaker.pdfinfoMaker_cfi          import *
from CMS3.NtupleMaker.genJetMaker_cfi           import *
from CMS3.NtupleMaker.puSummaryInfoMaker_cfi    import *

from CMS3.NtupleMaker.jetFlavorMaker_cfi        import *

#This should be in the config
from CMS3.NtupleMaker.dilepGenFilter_cfi import dilepGenFilter


cms2GENSequence     = cms.Sequence(

  genMaker * pdfinfoMaker * genJetSequence * CMS2FlavorHistorySequence * candToGenAssMaker * genPFJetMaker * genJetMaker * hypGenMaker * puSummaryInfoMaker * jetFlavorSequence

)
