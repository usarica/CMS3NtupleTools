import FWCore.ParameterSet.Config as cms

printTree = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(1)
)

myPartons = cms.EDProducer("PartonSelector",
    withLeptons = cms.bool(False)
)

flavourByRef = cms.EDProducer("JetPartonMatcher",
    jets                = cms.InputTag("ak5CaloJets"),
    coneSizeToAssociate = cms.double(0.3),
    partons             = cms.InputTag("myPartons")
)

flavourByVal = cms.EDProducer("JetFlavourIdentifier",
    srcByReference    = cms.InputTag("flavourByRef"),
    physicsDefinition = cms.bool(True)
)


jetFlavorMaker = cms.EDProducer("JetFlavorMaker",

  srcSelectedPartons = cms.InputTag("myPartons"),
  srcByReference     = cms.InputTag("flavourByRef"),
  srcByValue         = cms.InputTag("flavourByVal")

)

#jetFlavorSequence = printTree * myPartons * flavourByRef * flavourByVal * genJetFlavorMaker
jetFlavorSequence = myPartons * flavourByRef * flavourByVal * jetFlavorMaker

