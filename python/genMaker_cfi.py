import FWCore.ParameterSet.Config as cms

#genMaker = cms.EDFilter("GenMaker",
genMaker = cms.EDProducer("GenMaker",
	aliasPrefix = cms.untracked.string("genps"),
    ntupleOnlyStatus3     = cms.bool(False), # just get everything since Pythia8 has weird status codes
    ntupleDaughters       = cms.bool(True),
    genParticlesInputTag  = cms.InputTag("prunedGenParticles" ),
    genEvtInfoInputTag  = cms.InputTag("generator" ),
    ntuplePackedGenParticles    = cms.bool(False), # default is False                          
    packedGenParticlesInputTag  = cms.InputTag("packedGenParticles" ), # Assign Status "1111" to these to avoid duplication. Only save p4, ID, status
    #genRunInfoInputTag    = cms.InputTag("generator"),
    exclusiveCrossSection = cms.untracked.double(0.0),
    inclusiveCrossSection = cms.untracked.double(0.0),
    kfactor               = cms.untracked.double(1.0),

    #PID of LSP, or MET particles besides nu's (which are always included)
    #For all LM points (msugra), LSP = 1000022. For GMSB, LSP = 1000039
    vmetPIDs = cms.untracked.vint32(1000022)
)

