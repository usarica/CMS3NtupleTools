import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.SISConeJetParameters_cfi import *

SISCone5TrkJets =  cms.EDProducer("SISConeJetProducer"                                       ,
                                  SISConeJetParameters                                       ,
                                  alias          = cms.untracked.string('SISCone05TrackJets'),
                                  inputEtMin     = cms.double(0.5)                           ,
                                  inputEMin      = cms.double(0.)                            ,
                                  src            = cms.InputTag('subTrkColl')                ,
                                  jetType        = cms.untracked.string("BasicJet")          ,
                                  jetPtMin       = cms.double(0.)                            ,
                                  UE_Subtraction = cms.string('no')                          ,
                                  coneRadius     = cms.double(0.5)
)

trkJetMaker = cms.EDProducer("TrkJetMaker", 
                             trkJetsInputTag = cms.InputTag('SISCone5TrkJets')
)

trkjetmaker = cms.Sequence(SISCone5TrkJets * trkJetMaker)
