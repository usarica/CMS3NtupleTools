import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Type1MET.MuonMETValueMapProducer_cff import *
from JetMETCorrections.Type1MET.MetMuonCorrections_cff import *
from JetMETCorrections.Type1MET.MuonMETValueMapProducer_cff import *
from JetMETCorrections.Type1MET.MuonTCMETValueMapProducer_cff import *
from JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff import *
from JetMETCorrections.Type1MET.MetType1Corrections_cff import corMetType1Icone5
from CMS2.NtupleMaker.metMaker_cfi import *
from RecoMET.METProducers.TCMET_cfi import *
from CMS2.NtupleMaker.tcmetMaker_cfi import *

metMuonJESCorSC5CMS2 = corMetType1Icone5.clone()
metMuonJESCorSC5CMS2.inputUncorJetsLabel = "prunedUncorrectedCMS2Jets"
metMuonJESCorSC5CMS2.corrector = "L2L3JetCorrectorSC5Calo"
metMuonJESCorSC5CMS2.inputUncorMetLabel = "corMetGlobalMuons"


metCorSequence = cms.Sequence(muonMETValueMapProducer*corMetGlobalMuons*metMuonJESCorSC5CMS2*metMaker*muonTCMETValueMapProducer*tcMet*tcmetMaker)

