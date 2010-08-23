import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.TCMET_cfi import *
from RecoMET.METProducers.MuonTCMETValueMapProducer_cff import *
from CMS2.NtupleMaker.tcmetMaker_cfi import *

muonTCMETValueMapProducerNew = muonTCMETValueMapProducer.clone()
muonTCMETValueMapProducerNew.muonMinValidStaHits = cms.int32(0)

tcMetNew = tcMet.clone()
tcMetNew.tcmetDepValueMap = cms.InputTag("muonTCMETValueMapProducerNew", "muCorrData")
tcMetNew.usePFClusters = cms.bool(True)

tcmetMakerNew = tcmetMaker.clone()
tcmetMakerNew.aliasPrefix = cms.untracked.string("evtnew")
tcmetMakerNew.tcmet_tag_ = cms.InputTag("tcMetNew")
tcmetMakerNew.tcmet_vm_tag_ = cms.InputTag("muonTCMETValueMapProducerNew", "muCorrData")

tcmetSequence = cms.Sequence(muonTCMETValueMapProducerNew * tcMetNew * tcmetMakerNew)
