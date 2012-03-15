from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
offlinePrimaryVerticesDA = offlinePrimaryVertices.clone()

# What is commented is default in & after 520pre6


# 
#TkFilterParameters=cms.PSet(
#  algorithm                 = cms.string('filter'),
#  maxNormalizedChi2         = cms.double(20.0),     # loose cut on track chi**2
#  minPixelLayersWithHits    = cms.int32(2),         # two or more pixel layers
#  minSiliconLayersWithHits  = cms.int32(5),         # five or more tracker layers (includes pixels)
#  maxD0Significance         = cms.double(5.0),      # impact parameter significance
#  minPt                     = cms.double(0.0),      # direct pt cut
#  trackQuality              = cms.string("any")     # track quality not used
#)

#
TkClusParameters = cms.PSet(
    algorithm=cms.string('DA'),
    TkDAClusParameters = cms.PSet(
      coolingFactor = cms.double(0.6),  # slow annealing
      Tmin          = cms.double(4.0),  # freezeout temperature
      vertexSize    = cms.double(0.01), # ~ resolution / sqrt(Tmin)
      d0CutOff      = cms.double(3.),   # downweight high IP tracks 
      dzCutOff      = cms.double(4.)    # outlier rejection after freeze-out
    )
)

##
#vertexCollections = cms.VPSet(
#  [
#    cms.PSet(
#      label             = cms.string(""),
#      algorithm         = cms.string("AdaptiveVertexFitter"),
#      minNdof           = cms.double(0.0),
#      useBeamConstraint = cms.bool(False),
#      maxDistanceToBeam = cms.double(1.0)
#    )
#  ,
#    cms.PSet(label      = cms.string("WithBS"),
#      algorithm         = cms.string('AdaptiveVertexFitter'),
#      minNdof           = cms.double(2.0),
#      useBeamConstraint = cms.bool(True),
#      maxDistanceToBeam = cms.double(1.0)
#    )
#  ]
#)

#offlinePrimaryVerticesDA.verbose            = cms.untracked.bool(False)      
#offlinePrimaryVerticesDA.TrackLabel         = cms.InputTag("generalTracks")
#offlinePrimaryVerticesDA.beamSpotLabel      = cms.InputTag("offlineBeamSpot")
#offlinePrimaryVerticesDA.TkFilterParameters = TkFilterParameters
offlinePrimaryVerticesDA.TkClusParameters   = TkClusParameters
#offlinePrimaryVerticesDA.vertexCollections  = vertexCollections

# Not Clear why these need to be set in 52X...
offlinePrimaryVerticesDA.algorithm         = cms.string("AdaptiveVertexFitter")
offlinePrimaryVerticesDA.useBeamConstraint = cms.bool( False )
offlinePrimaryVerticesDA.PVSelParameters   = cms.PSet( maxDistanceToBeam = cms.double(1.0) )
offlinePrimaryVerticesDA.minNdof           = cms.double(0.0)

# These lines should be uncommented only for the gcc4X builds, with X>=6
#offlinePrimaryVerticesDA.TkClusParameters.algorithm = cms.string("DA_vect" )
#offlinePrimaryVerticesDA.TkClusParameters.TkDAClusParameters.use_vdt = cms.untracked.bool( True )


davertexreco = cms.Sequence(offlinePrimaryVerticesDA)
