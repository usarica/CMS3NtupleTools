import FWCore.ParameterSet.Config as cms


from PhysicsTools.HepMCCandAlgos.flavorHistoryProducer_cfi import *
from PhysicsTools.HepMCCandAlgos.flavorHistoryFilter_cfi import *
from CMS2.NtupleMaker.flavorHistoryMaker_cfi import *


#needs the gen jets thats produced by the genJetMaker
cFlavorHistoryProducer.matchedSrc = cms.InputTag("sisCone5StGenJets")
bFlavorHistoryProducer.matchedSrc = cms.InputTag("sisCone5StGenJets")


CMS2FlavorHistorySequence = cms.Sequence(bFlavorHistoryProducer * cFlavorHistoryProducer * flavorHistoryFilter * flavorHistoryMaker)
