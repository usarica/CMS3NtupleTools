import FWCore.ParameterSet.Config as cms

from CMS3.NtupleMaker.jetCollectionPruner_cfi import *

cms2JetSequence = cms.Sequence(prunedUncorrectedCMS2Jets)
