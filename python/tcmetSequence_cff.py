import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Type1MET.MuonMETValueMapProducer_cff import *
from JetMETCorrections.Type1MET.MuonTCMETValueMapProducer_cff import *
from RecoMET.METProducers.TCMET_cfi import *
from CMS2.NtupleMaker.tcmetMaker_cfi import *

tcmetMakerOld = tcmetMaker.clone()
tcmetMakerOld.aliasPrefix = cms.untracked.string("evt35X")

muonMETValueMapProducerNew = muonMETValueMapProducer.clone()
muonTCMETValueMapProducerNew = muonTCMETValueMapProducer.clone()
tcMetNew = tcMet.clone()
tcMetNew.muonDepValueMap = cms.InputTag("muonMETValueMapProducerNew", "muCorrData")
tcMetNew.tcmetDepValueMap = cms.InputTag("muonTCMETValueMapProducerNew", "muCorrData")

tcmetMaker.tcmet_tag_ = cms.InputTag("tcMetNew")
tcmetMaker.tcmet_vm_tag_ = cms.InputTag("muonTCMETValueMapProducerNew", "muCorrData")

tcmetSequence = cms.Sequence(tcmetMakerOld * muonMETValueMapProducerNew * muonTCMETValueMapProducerNew * tcMetNew * tcmetMaker)

