import FWCore.ParameterSet.Config as cms

goodMusForMETCorr = cms.EDFilter("GoodMusForMETCorrProducer",
                         src = cms.InputTag("muons"),
                         isGlobalMuon = cms.bool(True),
                         ptCut = cms.double(10.0),
                         etaCut = cms.double(2.5),
                         numValidHitsCut = cms.int32(5),
                         qoverpErrorCut = cms.double(0.5)
)

