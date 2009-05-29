import FWCore.ParameterSet.Config as cms

from PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi import *
#process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")

flavorHistoryWbbMaker = cms.EDFilter(
	"FlavorHistMaker",
	pathNumber = cms.uint32(1)
	)

flavorHistoryWbMaker = cms.EDFilter(
	"FlavorHistMaker",
	pathNumber = cms.uint32(2)
	)

flavorHistoryWccMaker = cms.EDFilter(
	"FlavorHistMaker",
	pathNumber = cms.uint32(3)
	)

flavorHistoryWcMaker = cms.EDFilter(
	"FlavorHistMaker",
	pathNumber = cms.uint32(4)
	)

flavorHistoryWbbgsMaker = cms.EDFilter(
	"FlavorHistMaker",
	pathNumber = cms.uint32(5)
	)

flavorHistoryWccgsMaker = cms.EDFilter(
	"FlavorHistMaker",
	pathNumber = cms.uint32(6)
	)

flavorHistoryWjetsMaker = cms.EDFilter(
	"FlavorHistMaker",
	pathNumber = cms.uint32(11)
	)

cms2Wbb = cms.Path(
	genJetParticles*sisCone5GenJets*
	bFlavorHistoryProducer*
	cFlavorHistoryProducer*
	wbbMEFlavorHistoryFilter*	
	flavorHistoryWbbMaker)

cms2Wb = cms.Path(
	genJetParticles*sisCone5GenJets*
	bFlavorHistoryProducer*
	cFlavorHistoryProducer*
	~wbbMEFlavorHistoryFilter*
	wbFEFlavorHistoryFilter*
	flavorHistoryWbMaker)

cms2Wcc = cms.Path(
    genJetParticles*sisCone5GenJets*
	    bFlavorHistoryProducer*
	    cFlavorHistoryProducer*
	    ~wbbMEFlavorHistoryFilter*
	    ~wbFEFlavorHistoryFilter*
	    wccMEFlavorHistoryFilter*
	flavorHistoryWccMaker)

cms2Wc = cms.Path(
    genJetParticles*sisCone5GenJets*
	    bFlavorHistoryProducer*
	    cFlavorHistoryProducer*
	    ~wbbMEFlavorHistoryFilter*
	    ~wbFEFlavorHistoryFilter*
	    ~wccMEFlavorHistoryFilter*
	    wcFEFlavorHistoryFilter*
	flavorHistoryWcMaker)

cms2Wbb_gs = cms.Path(
    genJetParticles*sisCone5GenJets*
	    bFlavorHistoryProducer*
	    cFlavorHistoryProducer*
	    cFlavorHistoryProducer*
	    ~wbbMEFlavorHistoryFilter*
	    ~wbFEFlavorHistoryFilter*
	    ~wccMEFlavorHistoryFilter*
	    ~wcFEFlavorHistoryFilter*
	    wbbGSFlavorHistoryFilter*
	flavorHistoryWbbgsMaker)

cms2Wcc_gs = cms.Path(
    genJetParticles*sisCone5GenJets*
	    bFlavorHistoryProducer*
	    cFlavorHistoryProducer*
	    ~wbbMEFlavorHistoryFilter*
	    ~wbFEFlavorHistoryFilter*
	    ~wccMEFlavorHistoryFilter*
	    ~wcFEFlavorHistoryFilter*
	    ~wbbGSFlavorHistoryFilter*
	    wccGSFlavorHistoryFilter*
	flavorHistoryWccgsMaker)

cms2Wjets = cms.Path(
    genJetParticles*sisCone5GenJets*
	    bFlavorHistoryProducer*
	    cFlavorHistoryProducer*
	    ~wbbMEFlavorHistoryFilter*
	    ~wbFEFlavorHistoryFilter*
	    ~wccMEFlavorHistoryFilter*
	    ~wcFEFlavorHistoryFilter*
	    ~wbbGSFlavorHistoryFilter*
	    ~wccGSFlavorHistoryFilter*
	    ~wbbMEComplimentFlavorHistoryFilter*
	    ~wccMEComplimentFlavorHistoryFilter*
	    ~wbbGSComplimentFlavorHistoryFilter*
	    ~wccGSComplimentFlavorHistoryFilter*
	flavorHistoryWjetsMaker)

