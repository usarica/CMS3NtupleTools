import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *

trkJetMaker = cms.EDProducer("TrkJetMaker", 
                             trkJetsInputTag = cms.InputTag('prunedUncorrectedCMS2Jets','trkjet'),                             
                             trkJetCorrectionL2L3 = cms.string("ak5TRKL2L3")
)


