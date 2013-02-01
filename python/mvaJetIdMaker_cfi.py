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
		JetIdParams = JetIdParams,
                label = cms.string("philv1")

)


mvaJetIdMakerFull53x = cms.EDProducer("MVAJetIdMaker",
        CorrJetNameData = cms.InputTag("ak5PFJetsL1FastL2L3Residual"),
        CorrJetNameMC   = cms.InputTag("ak5PFJetsL1FastL2L3"),
        #CorrJetName    = cms.InputTag("ak5PFJets"),                    # for JEC test
        JetName         = cms.InputTag("ak5PFJets"),
        VertexName      = cms.InputTag("offlinePrimaryVertices"),
        cutBased        = cms.bool(False),
        JetPtMin        = cms.double(0.),
        tmvaWeights     = cms.string("CMGTools/External/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml"),
        tmvaMethod      = cms.string("JetIDMVAHighPt"),
        version         = cms.int32(-1),
        tmvaVariables   = cms.vstring(
    "nvtx"     ,
    "dZ"       ,
    "beta"     ,
    "betaStar" ,
    "nCharged" ,
    "nNeutrals",
    "dR2Mean"  ,
    "ptD"      ,
    "frac01"   ,
    "frac02"   ,
    "frac03"   ,
    "frac04"   ,
    "frac05"   , 
    ),
                               
        tmvaSpectators = cms.vstring(
    "jetPt",
    "jetEta",
    "jetPhi"
    ),
       JetIdParams = full_5x_wp,
       label = cms.string("full53x")                                      
                                      
)


mvaJetIdMakerFull5x = cms.EDProducer("MVAJetIdMaker",
        CorrJetNameData = cms.InputTag("ak5PFJetsL1FastL2L3Residual"),
        CorrJetNameMC   = cms.InputTag("ak5PFJetsL1FastL2L3"),
        #CorrJetName    = cms.InputTag("ak5PFJets"),                    # for JEC test
        JetName         = cms.InputTag("ak5PFJets"),
        VertexName      = cms.InputTag("offlinePrimaryVertices"),
        cutBased        = cms.bool(False),
        JetPtMin        = cms.double(0.),
        tmvaWeights     = cms.string("CMGTools/External/data/TMVAClassification_5x_BDT_fullPlusRMS.weights.xml"),
        tmvaMethod      = cms.string("BDT_fullPlusRMS"),
        version         = cms.int32(-1),
        tmvaVariables   = cms.vstring(
      "frac01",
      "frac02",
      "frac03",
      "frac04",
      "frac05",
      "dR2Mean",
      "nvtx",
      "nNeutrals",
      "beta",
      "betaStar",
      "dZ",
      "nCharged",
    ),
                               
        tmvaSpectators = cms.vstring(
    "jetPt",
    "jetEta",
    ),
       JetIdParams = full_5x_wp,
       label = cms.string("full5x")
)
