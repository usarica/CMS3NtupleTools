import FWCore.ParameterSet.Config as cms

from PhysicsTools.HepMCCandAlgos.genParticleCandidatesFast_cfi import *
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *
from RecoJets.Configuration.GenJetParticles_cff import *
# from RecoJets.JetProducers.AntiKtJetParameters_cfi import *
from RecoJets.JetProducers.ak5GenJets_cfi import *
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


cms2antikt5GenJets = ak5GenJets.clone()
cms2antikt5GenJets.src = cms.InputTag("genParticlesAllStables")
cms2antikt5GenJets.jetPtMin = cms.double(0.)
cms2antikt5GenJets.alias = cms.untracked.string("CMS2ANTIKT5GenJet")


genJetSequence = cms.Sequence( genParticlesAllStables * cms2antikt5GenJets )

