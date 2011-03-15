
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi import *
# track selection, common for all producers here (doesn't have to be)
TkFilterParameters=cms.PSet(
    algorithm=cms.string('filter'),
    minPt = cms.double(0.0),                   # direct pt cut
    maxD0Significance = cms.double(5.0),       # impact parameter significance
    maxNormalizedChi2 = cms.double(20.0),       # loose cut on track chi**2
    minPixelLayersWithHits = cms.int32(2),     # two or more pixel layers
    minSiliconLayersWithHits = cms.int32(5),   # five or more tracker layers (includes pixels)
    trackQuality = cms.string("any")           # track quality not used
    )
# offlinePrimaryVerticesDA   deterministic annealing clustering + fit with beam constraint
offlinePrimaryVerticesDA.verbose = cms.untracked.bool(False)
offlinePrimaryVerticesDA.TkFilterParameters=TkFilterParameters
offlinePrimaryVerticesDA.TrackLabel = cms.InputTag("generalTracks")
offlinePrimaryVerticesDA.minNdof  = cms.double(0.0)
offlinePrimaryVerticesDA.PVSelParameters=cms.PSet(   maxDistanceToBeam = cms.double(1.0)   )
offlinePrimaryVerticesDA.TkClusParameters=cms.PSet(
    algorithm=cms.string('DA'),
    TkDAClusParameters = cms.PSet(
      coolingFactor = cms.double(0.6),  #  slow annealing
      Tmin = cms.double(4.0),           #  freezeout temperature
      vertexSize = cms.double(0.01)     #  ~ resolution / sqrt(Tmin)
    )
)
davertexreco = cms.Sequence(offlinePrimaryVerticesDA)
