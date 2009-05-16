import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Type1MET.MuonMETValueMapProducer_cff import *
from JetMETCorrections.Type1MET.MetMuonCorrections_cff import *
from JetMETCorrections.Type1MET.MuonMETValueMapProducer_cff import *
from CMS2.NtupleMaker.metMaker_cfi import *
from JetMETCorrections.Type1MET.MuonTCMETValueMapProducer_cff import *
from RecoMET.METProducers.TCMET_cfi import *
from CMS2.NtupleMaker.tcmetMaker_cfi import *

metCorSequence = cms.Sequence(muonMETValueMapProducer*corMetGlobalMuons*metMaker*muonTCMETValueMapProducer*tcMet*tcmetMaker)

