import FWCore.ParameterSet.Config as cms
from CMGTools.External.puJetIDAlgo_cff import *

mvaJetIdMaker = cms.EDProducer("MVAJetIdMaker", 
	#CorrJetName     = cms.InputTag("ak5PFJetsL1L2L3"), 	# for JEC test
	#JetName         = cms.InputTag("ak5PFJets"),			# for JEC test
	CorrJetName     = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
	JetName         = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
    VertexName      = cms.InputTag("offlinePrimaryVertices"),	
	tmvaWeights 	= cms.untracked.string("CMGTools/External/data/mva_JetID_v1.weights.xml"),
	tmvaMethod  	= cms.untracked.string("JetID"),
	version 		= cms.untracked.int32(-1),
	tmvaVariables 	= cms.untracked.vstring(
		"nvtx",
		"jetPt",
		"jetEta",
		"jetPhi",
		"dZ",
		"d0",
		"beta",
		"betaStar",
		"nCharged",
		"nNeutrals",
		"dRMean",
		"frac01",
		"frac02",
		"frac03",
		"frac04",
		"frac05",
		),
	JetIdParams = JetIdParams
)
