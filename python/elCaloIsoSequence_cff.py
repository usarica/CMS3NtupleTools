import FWCore.ParameterSet.Config as cms

from EgammaAnalysis.EgammaIsolationProducers.egammaBasicClusterMerger_cfi import *
from CMS2.NtupleMaker.elCaloIsoMaker_cfi import *
egammaBasicClusterMerger.src = cms.VInputTag(cms.InputTag('hybridSuperClusters',  'hybridBarrelBasicClusters'), 
                                             cms.InputTag('multi5x5BasicClusters','multi5x5EndcapBasicClusters')
)
elCaloIsoSequence = cms.Sequence(egammaBasicClusterMerger*elCaloIsoMaker)
