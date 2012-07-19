import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.pfmetMaker_cfi import *
from JetMETCorrections.Type1MET.pfMETCorrections_cff import *

#producePFMETCorrections.remove(kt6PFJets)
#producePFMETCorrections.remove(ak5PFJets)

CMS2pfMetSequence = cms.Sequence(producePFMETCorrections * pfmetMaker)
