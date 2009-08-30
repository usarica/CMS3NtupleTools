import FWCore.ParameterSet.Config as cms

genMaker = cms.EDFilter("GenMaker",
    ntupleOnlyStatus3     = cms.bool(True),
    ntupleDaughters       = cms.bool(True),
    genParticlesInputTag  = cms.InputTag("genParticles" ),
    genEventScaleInputTag = cms.InputTag("genEventScale"),
    exclusiveCrossSection = cms.untracked.double(0.0),
    inclusiveCrossSection = cms.untracked.double(0.0),
    kfactor               = cms.untracked.double(1.0),

    #PID of LSP, or MET particles besides nu's (which are always included)
    #For all LM points (msugra), LSP = 1000022. For GMSB, LSP = 1000039
    vmetPIDs = cms.untracked.vint32(1000022)
)

