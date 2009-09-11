import FWCore.ParameterSet.Config as cms

from PhysicsTools.HepMCCandAlgos.genParticleCandidatesFast_cfi import *
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *
from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.SISConeJetParameters_cfi import *
from CMS2.NtupleMaker.genJetMaker_cfi import *

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

sisCone5StGenJets = cms.EDProducer("SISConeJetProducer",
                                   SISConeJetParameters,                                
                                   alias          = cms.untracked.string('SISCone05StGenJets'),
                                   inputEtMin     = cms.double(0.),
                                   inputEMin      = cms.double(0.),
                                   jetPtMin       = cms.double(0.),
                                   UE_Subtraction = cms.string('no'), 
                                   src            = cms.InputTag("genParticlesAllStables"),
                                   verbose        = cms.untracked.bool(False),
                                   jetType        = cms.untracked.string('GenJet'),
                                   coneRadius     = cms.double(0.5)
)                                   


genJetSequence = cms.Sequence( genParticlesAllStables + sisCone5StGenJets + genJetMaker )

