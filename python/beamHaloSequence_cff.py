import FWCore.ParameterSet.Config as cms

from RecoMET.Configuration.RecoMET_BeamHaloId_cff import *
from CMS2.NtupleMaker.beamHaloMaker_cfi import *

cms2beamHaloSequence = cms.Sequence(BeamHaloId * beamHaloMaker)


