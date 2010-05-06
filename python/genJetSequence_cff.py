import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoMuNoNu
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

cms2antikt5GenJets = ak5GenJets.clone(src = cms.InputTag("genParticlesForJetsNoMuNoNu") )
genJetSequence = cms.Sequence( genParticlesForJetsNoMuNoNu * cms2antikt5GenJets )

