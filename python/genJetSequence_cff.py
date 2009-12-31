import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

genParticlesAllStables =  cms.EDProducer("InputGenJetsParticleSelector",
          src                      = cms.InputTag("genParticles"),
          partonicFinalState       = cms.bool(False),
          excludeResonances        = cms.bool(False),
          excludeFromResonancePids = cms.vuint32(),
          tausAsJets               = cms.bool(False),
          ignoreParticleIDs        = cms.vuint32(
            1000022, 2000012, 2000014,
            2000016, 1000039, 5000039,
            4000012, 9900012, 9900014,
            9900016, 39, 12, 13, 14, 16)
)


cms2antikt5GenJets = ak5GenJets.clone(src = cms.InputTag("genParticlesAllStables") )
genJetSequence = cms.Sequence( genParticlesAllStables * cms2antikt5GenJets )

