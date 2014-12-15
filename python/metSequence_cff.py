import FWCore.ParameterSet.Config as cms

#commented out for moving to updated JEC, sergo 02/11/13
#from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorIC5CaloJet

from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from RecoJets.JetProducers.CaloTowerSchemeBWithHO_cfi import towerMakerWithHO
from RecoMET.METProducers.CaloMET_cfi import metHO


from CMS3.NtupleMaker.metMaker_cfi import *
from JetMETCorrections.Type1MET.caloMETCorrections_cff import *
from JetMETCorrections.Type1MET.correctedMet_cff import *
metType1Sequence = cms.Sequence(produceCaloMETCorrections * caloType1CorrectedMet * metMaker)


#commented out for moving to updated JEC, sergo 02/11/13
#metMuonJESCorAK5CMS2                     = metJESCorIC5CaloJet.clone()
#metMuonJESCorAK5CMS2.inputUncorJetsLabel = "ak5CaloJets"
#metMuonJESCorAK5CMS2.corrector           = "ak5CaloL2L3"
#metMuonJESCorAK5CMS2.inputUncorMetLabel  = "corMetGlobalMuons"

from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *
cms2HBHENoiseFilterResultProducer = HBHENoiseFilterResultProducer.clone()
#Updated configuration based on https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1196.html
#### REMOVE THE FOLLOWING SINCE THEY ARE OUTDATED
#cms2HBHENoiseFilterResultProducer.minIsolatedNoiseSumE = cms.double(999999.)
#cms2HBHENoiseFilterResultProducer.minNumIsolatedNoiseChannels = cms.int32(999999)
#cms2HBHENoiseFilterResultProducer.minIsolatedNoiseSumEt = cms.double(999999.)

metCorSequence = cms.Sequence(metType1Sequence * cms2HBHENoiseFilterResultProducer)
#metCorSequence = cms.Sequence(metMuonJESCorAK5CMS2 * towerMakerWithHO * metHO * cms2HBHENoiseFilterResultProducer)
