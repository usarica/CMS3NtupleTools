import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.TCMET_cfi import *
from CMS2.NtupleMaker.tcmetMaker_cfi import *

tcMetWithPFclusters = tcMet.clone()
tcMetWithPFclusters.tcmetDepValueMap = cms.InputTag("muonTCMETValueMapProducer", "muCorrData")
tcMetWithPFclusters.usePFClusters = cms.bool(True)

tcmetSequence = cms.Sequence( tcmetMaker )
#tcmetSequence = cms.Sequence(tcMetWithPFclusters * tcmetMaker)
