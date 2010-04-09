import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Type1MET.MuonMETValueMapProducer_cff import *
from JetMETCorrections.Type1MET.MuonTCMETValueMapProducer_cff import *
from RecoMET.METProducers.TCMET_cfi import *
from CMS2.NtupleMaker.metMaker_cfi import *
from CMS2.NtupleMaker.tcmetMaker_cfi import *

muonMETValueMapProducerNew = muonMETValueMapProducer.clone()
muonTCMETValueMapProducerNew = muonTCMETValueMapProducer.clone()
tcMetNew = tcMet.clone()
tcMetNew.muonDepValueMap = cms.InputTag("muonMETValueMapProducerNew", "muCorrData")
tcMetNew.tcmetDepValueMap = cms.InputTag("muonTCMETValueMapProducerNew", "muCorrData")
metMakerNew = metMaker.clone()
metMakerNew.aliasPrefix = cms.untracked.string("new")
metMakerNew.muon_vm_tag_ = cms.InputTag("muonMETValueMapProducerNew", "muCorrData")
tcmetMakerNew = tcmetMaker.clone()
tcmetMakerNew.aliasPrefix = cms.untracked.string("new")
tcmetMakerNew.tcmet_tag_ = cms.InputTag("tcMetNew")
tcmetMakerNew.tcmet_vm_tag_ = cms.InputTag("muonTCMETValueMapProducerNew", "muCorrData")

tcmetSequence = cms.Sequence(muonMETValueMapProducerNew * muonTCMETValueMapProducerNew * tcMetNew * metMakerNew * tcmetMakerNew)

