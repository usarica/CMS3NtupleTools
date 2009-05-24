
import FWCore.ParameterSet.Config as cms

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

localTrkColl = cms.EDProducer("TrkMuonFilter", 
                TrackProducerTag = cms.InputTag("generalTracks"),
                MuonTag = cms.InputTag("muons"),
                subMuon = cms.bool(True),
                muIsoFrac = cms.double(0.92),
                muChi2N  = cms.double(5),
                muMinPt = cms.double(0),
                muMaxEta = cms.double(2.5)
)


mTrkColl = cms.EDProducer("TrkMuonFilter", 
                TrackProducerTag = cms.InputTag("generalTracks"),
                MuonTag = cms.InputTag("muons"),
                subMuon = cms.bool(False),
                muIsoFrac = cms.double(0.92),
                muChi2N  = cms.double(5),
                muMinPt = cms.double(7),
                muMaxEta = cms.double(2.5)
)

subTrkColl = cms.EDProducer("ConcreteChargedCandidateProducer", 
        src =  cms.InputTag("localTrkColl"),
        particleType = cms.string('pi+')
)


trkColl = cms.EDProducer("ConcreteChargedCandidateProducer", 
        src =  cms.InputTag("mTrkColl"),
        particleType = cms.string('pi+')
)

trkmuonfilter = cms.Sequence(localTrkColl+subTrkColl)

