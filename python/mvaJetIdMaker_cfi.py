import FWCore.ParameterSet.Config as cms
from CMGTools.External.puJetIDAlgo_cff import *

mvaJetIdMaker = cms.EDProducer("MVAJetIdMaker", 
	CorrJetNameData = cms.InputTag("ak5PFJetsL1FastL2L3Residual"), 
	CorrJetNameMC   = cms.InputTag("ak5PFJetsL1FastL2L3"), 
	#CorrJetName    = cms.InputTag("ak5PFJets"),			# for JEC test
	JetName         = cms.InputTag("ak5PFJets"),		
	VertexName      = cms.InputTag("offlinePrimaryVertices"),	
        cutBased        = cms.bool(False),
        JetPtMin        = cms.double(0.),
	tmvaWeights   	= cms.string("CMGTools/External/data/mva_JetID_v1.weights.xml"),
	tmvaMethod    	= cms.string("JetID"),
	version       	= cms.int32(-1),
	tmvaVariables 	= cms.vstring(
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
		tmvaSpectators = cms.vstring(),
		JetIdParams = JetIdParams
)
