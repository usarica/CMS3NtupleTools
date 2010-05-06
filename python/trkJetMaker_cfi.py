import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *

trkJetMaker = cms.EDProducer("TrkJetMaker", 
                             trkJetsInputTag = cms.InputTag('prunedUncorrectedCMS2TrackJets'),
                             trkJetPtCut     = cms.double(3.0),
                             trkJetCorrectionL2L3 = cms.string("ak5TrackL2L3")
)


