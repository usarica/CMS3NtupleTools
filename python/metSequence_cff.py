import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Type1MET.MetMuonCorrections_cff import *
from CMS2.NtupleMaker.metMaker_cfi import *

metCorSequence = cms.Sequence(MetMuonCorrections*metMaker)

