import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoMuNoNu
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

cms2antikt5GenJets = ak5GenJets.clone(src = cms.InputTag("genParticlesForJetsNoMuNoNu") )
cms2antikt5PFGenJets = ak5GenJets.clone(src = cms.InputTag("genParticlesForJetsNoNu") )

genJetSequence = cms.Sequence( genParticlesForJetsNoMuNoNu * cms2antikt5GenJets * genParticlesForJetsNoNu * cms2antikt5PFGenJets)


