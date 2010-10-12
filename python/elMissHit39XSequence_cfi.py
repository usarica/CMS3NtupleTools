import FWCore.ParameterSet.Config as cms


expectedHitsEle = cms.EDProducer("ExpectedHitsComputer",
                                 inputColl       = cms.InputTag("gsfElectrons"),
                                 useGsfTrack     = cms.bool(True),
                                 objectSelection = cms.string(""),
                                 propagator         = cms.string('PropagatorWithMaterialOpposite'),
                                 navigationSchool   = cms.string('SimpleNavigationSchool'),
                                 measurementTracker = cms.string(''),
                                 )


elMissHit39XMaker = cms.EDProducer("ElectronMissHit39XMaker",
                                   aliasPrefix = cms.untracked.string("els"),
                                   electronsInputTag = cms.InputTag("gsfElectrons"),
                                   electronMissHit39XTag = cms.InputTag("expectedHitsEle"),
                                   )

elMissHit39XSequence = cms.Sequence(expectedHitsEle * elMissHit39XMaker)

