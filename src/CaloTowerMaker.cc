// -*- C++ -*-
//
// Package:    CaloTowerMaker
// Class:      CaloTowerMaker
// 
/**\class CaloTowerMaker CaloTowerMaker.cc CMS2/NtupleMaker/src/CaloTowerMaker.cc

Description: <produce TaS collection of CaloTowers>

Implementation:
<Currently a bare copy of SCMaker>
 */
//
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "CMS2/NtupleMaker/interface/CaloTowerMaker.h"
//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalCaloFlagLabels.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;

using namespace std;

//
// class decleration
//

//
// constructors and destructor
//
CaloTowerMaker::CaloTowerMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

	// number of towers in the event
	produces<unsigned int>("evtn"+branchprefix).setBranchAlias("evt_n"+branchprefix);

	produces<std::vector<float> >(branchprefix+"eta").setBranchAlias(aliasprefix_+"_eta");
	produces<std::vector<float> >(branchprefix+"phi").setBranchAlias(aliasprefix_+"_phi");
	produces<std::vector<uint32_t> >(branchprefix+"detid").setBranchAlias(aliasprefix_+"_detid");

	// energy contributions from different detectors
	// energy in HO ("outerEnergy")is not included in "hadEnergy"
	//   double emEnergy() const { return emE_ ; }
	produces<std::vector<float> >(branchprefix+"emEnergy").setBranchAlias(aliasprefix_+"_emEnergy");
	//   double hadEnergy() const { return hadE_ ; }
	produces<std::vector<float> >(branchprefix+"hadEnergy").setBranchAlias(aliasprefix_+"_hadEnergy");
	//   double outerEnergy() const { return (id_.ietaAbs()<16)? outerE_ : 0.0; }
	produces<std::vector<float> >(branchprefix+"outerEnergy").setBranchAlias(aliasprefix_+"_outerEnergy");

	// transverse energies wrt to vtx (0,0,0)
	//   double emEt() const { return emE_ * sin( theta() ); }
	produces<std::vector<float> >(branchprefix+"emEt").setBranchAlias(aliasprefix_+"_emEt");
	//   double hadEt() const { return hadE_ * sin( theta() ); }
	produces<std::vector<float> >(branchprefix+"hadEt").setBranchAlias(aliasprefix_+"_hadEt");
	//   double outerEt() const { return (id_.ietaAbs()<16)? outerE_ * sin( theta() ) : 0.0; }
	produces<std::vector<float> >(branchprefix+"outerEt").setBranchAlias(aliasprefix_+"_outerEt");

	// recalculated wrt vertex provided as 3D point
	//   math::PtEtaPhiMLorentzVector p4(Point v) const;
	//   double p (Point v) const { return p4(v).P(); }
	produces<std::vector<float> >(branchprefix+"pcorr").setBranchAlias(aliasprefix_+"_pcorr");
	//   double et(Point v) const { return p4(v).Et(); }
	produces<std::vector<float> >(branchprefix+"etcorr").setBranchAlias(aliasprefix_+"_etcorr");
	//   double emEt(Point v)  const { return  emE_ * sin(p4(v).theta()); }
	produces<std::vector<float> >(branchprefix+"emEtcorr").setBranchAlias(aliasprefix_+"_emEtcorr");
	//   double hadEt(Point v) const { return  hadE_ * sin(p4(v).theta()); }
	produces<std::vector<float> >(branchprefix+"hadEtcorr").setBranchAlias(aliasprefix_+"_hadEtcorr");
	//   double outerEt(Point v) const { return (id_.ietaAbs()<16)? outerE_ * sin(p4(v).theta()) : 0.0; }
	produces<std::vector<float> >(branchprefix+"outerEtcorr").setBranchAlias(aliasprefix_+"_outerEtcorr");
	produces<std::vector<float> >(branchprefix+"etacorr").setBranchAlias(aliasprefix_+"_etacorr");
	produces<std::vector<float> >(branchprefix+"phicorr").setBranchAlias(aliasprefix_+"_phicorr");


	//    // the reference poins in ECAL and HCAL for direction determination
	//    // algorithm and parameters for selecting these points are set in the CaloTowersCreator
	//    const GlobalPoint& emPosition()  const { return emPosition_ ; }
	//    const GlobalPoint& hadPosition() const { return hadPosition_ ; }

	// time (ns) in ECAL/HCAL components of the tower based on weigted sum of the times in the contributing RecHits
	//   float ecalTime() const { return float(ecalTime_) * 0.01; }
	produces<std::vector<float> >(branchprefix+"ecalTime").setBranchAlias(aliasprefix_+"_ecalTime");
	//   float hcalTime() const { return float(hcalTime_) * 0.01; }
	produces<std::vector<float> >(branchprefix+"hcalTime").setBranchAlias(aliasprefix_+"_hcalTime");
	// hcal time is basically junk, but the hit time is not. For HBHE and HF, store it:
	produces<vector<vector<float> > >(branchprefix+"hcalHitTime").setBranchAlias(aliasprefix_+"_hcalHitTime");
	// HBHE and HF depth: 2 = short, 1 = long
	produces<vector<vector<int  > > >(branchprefix+"hcalHitDepth").setBranchAlias(aliasprefix_+"_hcalHitDepth");
	//see RecoLocalCalo/HcalRecAlgos/interface/HcalCaloFlagLabels.h for below
	produces<vector<vector<int> > >(branchprefix+"hcalHitFlag").setBranchAlias(aliasprefix_+"_hcalHitFlag");

	// methods to retrieve status information from the CaloTower:
	// number of bad/recovered/problematic cells in the tower
	// separately for ECAL and HCAL
	//  uint numBadEcalCells() const { return (twrStatusWord_ & 0x1F); }
	produces<std::vector<unsigned int> >(branchprefix+"numBadEcalCells").setBranchAlias(aliasprefix_+"_numBadEcalCells");
	//  uint numRecoveredEcalCells() const { return ((twrStatusWord_ >> 5) & 0x1F); }
	produces<std::vector<unsigned int> >(branchprefix+"numRecoveredEcalCells").setBranchAlias(aliasprefix_+"_numRecoveredEcalCells");
	//  uint numProblematicEcalCells() const { return ((twrStatusWord_ >> 10) & 0x1F); }
	produces<std::vector<unsigned int> >(branchprefix+"numProblematicEcalCells").setBranchAlias(aliasprefix_+"_numProblematicEcalCells");

	//  uint numBadHcalCells() const { return ( (twrStatusWord_ >> 15)& 0x7); }
	produces<std::vector<unsigned int> >(branchprefix+"numBadHcalCells").setBranchAlias(aliasprefix_+"_numBadHcalCells");
	//  uint numRecoveredHcalCells() const { return ((twrStatusWord_ >> 18) & 0x7); }
	produces<std::vector<unsigned int> >(branchprefix+"numRecoveredHcalCells").setBranchAlias(aliasprefix_+"_numRecoveredHcalCells");
	//  uint numProblematicHcalCells() const { return ((twrStatusWord_ >> 21) & 0x7); }
	produces<std::vector<unsigned int> >(branchprefix+"numProblematicHcalCells").setBranchAlias(aliasprefix_+"_numProblematicHcalCells");

	// vector of 10 samples per max crystal in the calo tower
	produces<std::vector<std::vector<int> > >	(branchprefix+"emMaxEcalMGPASampleADC").setBranchAlias(aliasprefix_+"_emMaxEcalMGPASampleADC");
	//number of crystals 
	produces<std::vector<int> >		(branchprefix+"numCrystals").setBranchAlias(aliasprefix_+"_numCrystals");

	//The following set of branches are vectors over ALL ECAL HITS WITH ET > 5 (5 is configurable in attentionEtThresh_)

	// chi2 prob -- NO LONGER ACCESSIBLE IN 3_5_5. Replace this branch with Chi2
	//produces<std::vector<float> >		(branchprefix+"emThreshChi2Prob").setBranchAlias(aliasprefix_+"_emThreshChi2Prob");
	produces<std::vector<std::vector<float> > >		(branchprefix+"emThreshChi2").setBranchAlias(aliasprefix_+"_emThreshChi2");
	// time of crystals with em et > 5
	produces<std::vector<std::vector<float> > >		(branchprefix+"emThreshTime").setBranchAlias(aliasprefix_+"_emThreshTime");
	// the recoflag of these crystals
	produces<std::vector<std::vector<int> > >		(branchprefix+"emThreshRecoFlag").setBranchAlias(aliasprefix_+"_emThreshRecoFlag");
	// the severity level of the hit 	(see RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h)
	produces<std::vector<std::vector<int> > >		(branchprefix+"emThreshSevLvl").setBranchAlias(aliasprefix_+"_emThreshSevLvl");
	// the energy of the hit
	produces<std::vector<std::vector<float> > >		(branchprefix+"emThresh").setBranchAlias(aliasprefix_+"_emThresh");
	// the eta of the hit
	produces<std::vector<std::vector<float> > >		(branchprefix+"emThreshEta").setBranchAlias(aliasprefix_+"_emThreshEta");
	// the energy in 3x3 crystals centred on the max energy crystal
	produces<std::vector<std::vector<float> > >		(branchprefix+"em3x3").setBranchAlias(aliasprefix_+"_em3x3");
	// as above for 5x5 crystals
	produces<std::vector<std::vector<float> > >		(branchprefix+"em5x5").setBranchAlias(aliasprefix_+"_em5x5");
	// swiss cross
	produces<std::vector<std::vector<float> > >		(branchprefix+"emSwiss").setBranchAlias(aliasprefix_+"_emSwiss");

	// These two branches are vectors OVER SPIKES NOT OVER TOWERS (empty if no spikes in the event)
	// spikes are identified by R4 < spikeR4Thresh_ and Et > spikeEtThresh_ and spikeEtaMax_ (from config)
	produces<std::vector<float> >(branchprefix+"spikeEt").setBranchAlias(aliasprefix_+"_spikeEt");
	produces<std::vector<float> >(branchprefix+"spikeR4").setBranchAlias(aliasprefix_+"_spikeR4");

	// add superclusters to the ntuple if they have ET > scEtMin_
	//   scEtMin_ = iConfig.getParameter<double>("scEtMin");

	// input Tags
	primaryVertexInputTag_ = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
	caloTowersInputTag_ = iConfig.getParameter<edm::InputTag>("caloTowersInputTag");
	ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
	ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");
	ecalDigiProducerEE_     = iConfig.getParameter<edm::InputTag>("ecalDigiProducerEE");
	ecalDigiProducerEB_     = iConfig.getParameter<edm::InputTag>("ecalDigiProducerEB");
	hbheRecHitsInputTag_    = iConfig.getParameter<edm::InputTag>("hbheRecHitsInputTag");
	hfRecHitsInputTag_      = iConfig.getParameter<edm::InputTag>("hfRecHitsInputTag");

	threshHcal_     = iConfig.getParameter<double>("threshHcal");
	threshEt_       = iConfig.getParameter<double>("threshEt");
	spikeEtThresh_  = iConfig.getParameter<double>("spikeEtThresh");
	spikeR4Thresh_  = iConfig.getParameter<double>("spikeR4Thresh");
	spikeEtaMax_    = iConfig.getParameter<double>("spikeEtaMax");

	// initialise this
	//  cachedCaloGeometryID_ = 0;

}

void CaloTowerMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	// calo topology
	edm::ESHandle<CaloTopology> pTopology;
	iSetup.get<CaloTopologyRecord>().get(pTopology);
	topology_ = pTopology.product();

	//   // get the calo geometry
	//   if (cachedCaloGeometryID_ != iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
	//     cachedCaloGeometryID_ = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
	//     iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
	//   }

	// get the primary vertices
	edm::Handle<reco::VertexCollection> vertexHandle;
	try {
		iEvent.getByLabel(primaryVertexInputTag_, vertexHandle);
	}
	catch ( cms::Exception& ex ) {
		edm::LogError("CaloTowerMakerError") << "Error! can't get the primary vertex";
	}
	const reco::VertexCollection *vertexCollection = vertexHandle.product();
	Point pv(0.0, 0.0, 0.0);
	if (vertexCollection->size() > 0) {
		pv = vertexCollection->at(0).position();
	}

	edm::Handle<CaloTowerCollection> calotower;
	iEvent.getByLabel(caloTowersInputTag_,calotower);

	if(!calotower.isValid()) {
		edm::LogError("CaloTowerMakerError") << "Error! Can't get calotowers!" << std::endl;
		//    return ;
	}

	// get hits
	edm::Handle<EcalRecHitCollection> rhcHandleEE;
	iEvent.getByLabel(ecalRecHitsInputTag_EE_, rhcHandleEE);
	const EcalRecHitCollection *recHitsEE = rhcHandleEE.product();

	edm::Handle<EcalRecHitCollection> rhcHandleEB;
	iEvent.getByLabel(ecalRecHitsInputTag_EB_, rhcHandleEB);
	const EcalRecHitCollection *recHitsEB = rhcHandleEB.product();

	// ECAL Digis
	edm::Handle<EBDigiCollection> ebDigiHandle;
	iEvent.getByLabel(ecalDigiProducerEB_, ebDigiHandle);
	edm::Handle<EEDigiCollection> eeDigiHandle;
	iEvent.getByLabel(ecalDigiProducerEE_, eeDigiHandle);

	const EBDigiCollection *ebDigis = 0;
	const EEDigiCollection *eeDigis = 0;
	digi_ = false;
	if (ebDigiHandle.isValid() && eeDigiHandle.isValid()) {
		digi_ = true;
		ebDigis = ebDigiHandle.product();
		eeDigis = eeDigiHandle.product();
	}

	//ecal channel status
	edm::ESHandle<EcalChannelStatus> chStatus;
	iSetup.get<EcalChannelStatusRcd>().get(chStatus);
	theEcalChStatus_ = chStatus.product();

	//hf/hbhe stuffs
	edm::Handle<HFRecHitCollection> hf_rechit;
    iEvent.getByLabel(hfRecHitsInputTag_, hf_rechit);
	edm::Handle<HBHERecHitCollection> hbhe_rechit;
    iEvent.getByLabel(hbheRecHitsInputTag_, hbhe_rechit);


	// ecal cluster shape variables
	// do not use the lazy tools because need to get the hits anyway
	EcalClusterTools clusterTools;

	//   // get hoe variable
	//   HoECalculator hoeCalc(caloGeometry_);

	std::auto_ptr<unsigned int> evt_ntwrs (new unsigned int);

	std::auto_ptr<std::vector<float> > vector_twrs_eta (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_phi (new std::vector<float>);
	std::auto_ptr<std::vector<uint32_t> > vector_twrs_detid (new std::vector<uint32_t>);

	std::auto_ptr<std::vector<float> > vector_twrs_emEnergy (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_hadEnergy (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_outerEnergy (new std::vector<float>);

	std::auto_ptr<std::vector<float> > vector_twrs_emEt (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_hadEt (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_outerEt (new std::vector<float>);

	std::auto_ptr<std::vector<float> > vector_twrs_pcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_etcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_emEtcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_hadEtcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_outerEtcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_etacorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_phicorr (new std::vector<float>);

	std::auto_ptr<std::vector<float> > vector_twrs_ecalTime (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_hcalTime (new std::vector<float>);
	std::auto_ptr<vector<vector<float> > > vector_twrs_hcalHitTime	(new vector<vector<float> >);
	std::auto_ptr<vector<vector<int  > > > vector_twrs_hcalHitDepth	(new vector<vector<int> >);
	std::auto_ptr<vector<vector<int  > > > vector_twrs_hcalHitFlag (new vector<vector<int> >);

	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numBadEcalCells (new std::vector<unsigned int>);
	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numRecoveredEcalCells (new std::vector<unsigned int>);
	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numProblematicEcalCells (new std::vector<unsigned int>);

	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numBadHcalCells (new std::vector<unsigned int>);
	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numRecoveredHcalCells (new std::vector<unsigned int>);
	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numProblematicHcalCells (new std::vector<unsigned int>);

	std::auto_ptr<std::vector<std::vector<int> > > vector_twrs_emMaxEcalMGPASampleADC (new std::vector<std::vector<int> >);
	std::auto_ptr<std::vector<int>   > vector_twrs_numCrystals(new std::vector<int>);

	auto_ptr<vector<vector<float> > > vector_twrs_emThreshChi2 			(new vector<vector<float> >);
	auto_ptr<vector<vector<float> > > vector_twrs_emThreshChi2Prob 		(new vector<vector<float> >);
	auto_ptr<vector<vector<float> > > vector_twrs_emThreshTime      	(new vector<vector<float> >);
	auto_ptr<vector<vector<int  > > > vector_twrs_emThreshRecoFlag      (new vector<vector<int  > >);
	auto_ptr<vector<vector<int  > > > vector_twrs_emThreshSevLvl      	(new vector<vector<int  > >);
	auto_ptr<vector<vector<float> > > vector_twrs_emThresh      		(new vector<vector<float> >);
	auto_ptr<vector<vector<float> > > vector_twrs_emThreshEta      		(new vector<vector<float> >);
	auto_ptr<vector<vector<float> > > vector_twrs_em3x3      			(new vector<vector<float> >);
	auto_ptr<vector<vector<float> > > vector_twrs_em5x5      			(new vector<vector<float> >);
	auto_ptr<vector<vector<float> > > vector_twrs_emSwiss				(new vector<vector<float> >);	

	//vectors over spikes, not towers
	std::auto_ptr<std::vector<float>   > vector_twrs_spikeEt  (new std::vector<float>);
	std::auto_ptr<std::vector<float>   > vector_twrs_spikeR4  (new std::vector<float>);

	*evt_ntwrs = 0;

	for (CaloTowerCollection::const_iterator j = calotower->begin();j !=calotower->end(); j++) {

		*evt_ntwrs +=1;

		//     std::cout << *j << std::endl;
		//     std::cout << "ENERGY HAD " << j->hadEnergy()<< " ENERGY EM " <<j->emEnergy() 
		//          << " ETA " <<j->eta() << " PHI " <<j->phi() << std::endl;

		vector_twrs_eta->push_back(j->eta());
		vector_twrs_phi->push_back(j->phi());
		vector_twrs_detid->push_back(j->id().rawId());

		vector_twrs_emEnergy->push_back(j->emEnergy());
		vector_twrs_hadEnergy->push_back(j->hadEnergy());
		vector_twrs_outerEnergy->push_back(j->outerEnergy());

		vector_twrs_emEt->push_back(j->emEt());
		vector_twrs_hadEt->push_back(j->hadEt());
		vector_twrs_outerEt->push_back(j->outerEt());

		vector_twrs_pcorr->push_back(j->p(pv));
		vector_twrs_etcorr->push_back(j->et(pv));
		vector_twrs_emEtcorr->push_back(j->emEt(pv));
		vector_twrs_hadEtcorr->push_back(j->hadEt(pv));
		vector_twrs_outerEtcorr->push_back(j->outerEt(pv));
		vector_twrs_etacorr->push_back(j->p4(pv).eta());
		vector_twrs_phicorr->push_back(j->p4(pv).phi());

		vector_twrs_ecalTime->push_back(j->ecalTime());
		vector_twrs_hcalTime->push_back(j->hcalTime());

		vector_twrs_numBadEcalCells->push_back(j->numBadEcalCells());
		vector_twrs_numRecoveredEcalCells->push_back(j->numRecoveredEcalCells());
		vector_twrs_numProblematicEcalCells->push_back(j->numProblematicEcalCells());

		vector_twrs_numBadHcalCells->push_back(j->numBadHcalCells());
		vector_twrs_numRecoveredHcalCells->push_back(j->numRecoveredHcalCells());
		vector_twrs_numProblematicHcalCells->push_back(j->numProblematicHcalCells());

		// find the detids of the crystals in this calo tower
		const std::vector<DetId> &towerDetIds = j->constituents();

		// get variables for highest em energy crystal in tower
		float emE = 0.0;
		float emMax = 0.0;
		DetId emMaxId(0);
		vector<int> ecalMGPASampleADC;
		//below are for Thresh branches--for et > threshEt
		vector<float> chi2;
		vector<float> chi2Prob;
		vector<float> emTime;
		vector<int  > recoFlag;
		vector<int  > sevlvl;
		vector<float> emThresh;
		vector<float> emThreshEta;
		vector<float> em3x3;
		vector<float> em5x5;
		vector<float> emSwiss;
		vector<DetId> emId;
		vector<float> hcalTime;
		vector<int  > hcalDepth;
		vector<int  > hcalFlag;
		
		// loop on detids in the tower
		for (size_t i = 0; i < towerDetIds.size(); ++i) {
		  // find the energy of this detId if it is in the ecal
		  // also find spikes here for every hit in barrel towers, and fill numSpikes.
		  if (towerDetIds[i].det() == DetId::Ecal && towerDetIds[i].subdetId() == EcalEndcap) { //no spikes in endcap
			emE = clusterTools.recHitEnergy(towerDetIds[i], recHitsEE);
			//for endcap, use eta of tower, not hit. Sorry, this is a bit of a fudge, but i can't easily get eta of an endcap hit
			if( emE/cosh(j->eta()) > threshEt_ )
			  emId.push_back( towerDetIds[i] );
		  }
		  else if (towerDetIds[i].det() == DetId::Ecal && towerDetIds[i].subdetId() == EcalBarrel) {
			emE = clusterTools.recHitEnergy(towerDetIds[i], recHitsEB);
			//this is from RecoLocalCalo/EcalRecAlgos/src/EcalSeverityLevelAlgo.cc, or something
			float approxEta = EBDetId::approxEta( towerDetIds[i] );
			if( emE/approxEta > threshEt_ )
			  emId.push_back( towerDetIds[i] );
			//spike
			if( emE/cosh(approxEta) > spikeEtThresh_ ) { //check et cut first to skip SwissCross for run time speedup
			  reco::BasicCluster dummyCluster;
			  float s4 = SwissCross(dummyCluster, recHitsEB, towerDetIds[i]) - emE; //swiss cross still contains center hit
			  if( s4/emE < spikeR4Thresh_ && //identify spikes with r4
				  fabs(approxEta) < spikeEtaMax_ ) { //exclude transition
				vector_twrs_spikeEt->push_back( emE/cosh(approxEta) );
				vector_twrs_spikeR4->push_back( s4/emE );
			  }
			}
		  }

		  // compare with previous highest energy
		  if (emE > emMax) {
			emMax = emE; 
			emMaxId = towerDetIds[i];
		  }

		  //hf flags
		  //note that check of subdetId is missing (not needed?), and thresh is on em+had bc want both hits of hf
		  if( towerDetIds[i].det() == DetId::Hcal && j->emEt() + j->hadEt() > threshHcal_ ) {
			if( hf_rechit.isValid() && hbhe_rechit.isValid() ) {
			  HFRecHitCollection::const_iterator hfit = hf_rechit->find(towerDetIds[i]);
			  HBHERecHitCollection::const_iterator hbheit = hbhe_rechit->find(towerDetIds[i]);
			  //hcal calo flag labels
			  if( hfit != hf_rechit->end() ) {
				hcalFlag.push_back( hfit->flags() );
				hcalTime.push_back( hfit->time() );
				hcalDepth.push_back( hfit->id().depth() );
				//cout << j-calotower->begin() << "  " << hfit-hf_rechit->begin()     << "  " << j->eta() << "  " << hfit->id().depth()   << "  " << j->hcalTime() << "  " << hfit->time() << "  " << hfit->flags() << endl;
			  }
			  else if( hbheit != hbhe_rechit->end() ) {
				hcalFlag.push_back( hbheit->flags() );
				hcalTime.push_back(  hbheit->time() );
				hcalDepth.push_back( hbheit->id().depth() );
				//cout << j-calotower->begin() << "  " << hbheit-hbhe_rechit->begin() << "  " << j->eta() << "  " << hbheit->id().depth() << "  " << j->hcalTime() << "  " << hbheit->time() << "  " << hbheit->flags() << endl;
			  }
			}
			else {
			  cout << "One of either hbhe or hf rechit collections are bad. Check calotowermaker cfg" << endl;
			  exit(1);
			}
		  }

		} //end loop on towerDetIds

		//if none > threshEt_, store in Thresh branches the max anyway (so long as there is a max)
		if( emId.size() == 0 && emMaxId != DetId(0))
		  emId.push_back(emMaxId);

		// find the relevant quantities for the identified crystals
		for( unsigned int i=0; i<emId.size(); i++) {
		  if (emId[i] != DetId(0)) {
			reco::BasicCluster dummyCluster;
			if (emId[i].subdetId() == EcalEndcap) {
			  //chi2Prob.push_back( recHitChi2Prob(emId[i], recHitsEE) ); //cannot use in 3_5_5
			  chi2.push_back( recHitChi2(emId[i], recHitsEE) ); //replace above with this
			  emTime.push_back( recHitTime(emId[i], recHitsEE) );
			  sevlvl.push_back( recHitSeverityLevel(emId[i], recHitsEE) );
			  recoFlag.push_back( recHitFlag(emId[i], recHitsEE) );
			  emThresh.push_back( clusterTools.recHitEnergy(emId[i], recHitsEE) );
			  emThreshEta.push_back( j->eta() ); //again, sorry, this is just the tower eta, not hit eta
			  em3x3.push_back( clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, emId[i], -1, 1, -1, 1) );
			  em5x5.push_back( clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, emId[i], -2, 2, -2, 2) );
			  emSwiss.push_back( SwissCross(dummyCluster, recHitsEE, emId[i]) ); // make the swiss cross
			
			  if( emId[i] == emMaxId ) //adc just for max
				recHitSamples(emMaxId, eeDigis, ecalMGPASampleADC);
			}
			else if (emId[i].subdetId() == EcalBarrel) { 
			  //chi2Prob.push_back( recHitChi2Prob(emId[i], recHitsEB) ); //cannot use in 3_5_5
			  chi2.push_back( recHitChi2(emId[i], recHitsEB) ); //replace above with this
			  emTime.push_back( recHitTime(emId[i], recHitsEB) );
			  sevlvl.push_back( recHitSeverityLevel(emId[i], recHitsEB) );
			  recoFlag.push_back( recHitFlag(emId[i], recHitsEB) );
			  emThresh.push_back( clusterTools.recHitEnergy(emId[i], recHitsEB) );
			  emThreshEta.push_back( EBDetId::approxEta(emId[i]) );
			  em3x3.push_back( clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, emId[i], -1, 1, -1, 1) );
			  em5x5.push_back( clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, emId[i], -2, 2, -2, 2) );
			  emSwiss.push_back( SwissCross(dummyCluster, recHitsEB, emId[i]) ); // make the swiss cross

			  if( emId[i] == emMaxId ) //adc just for max
				recHitSamples(emMaxId, ebDigis, ecalMGPASampleADC);
			}
		  }
		}

		vector_twrs_emMaxEcalMGPASampleADC->push_back(ecalMGPASampleADC);
		vector_twrs_numCrystals->push_back(j->numCrystals());

		vector_twrs_emThreshChi2->push_back(chi2);
		vector_twrs_emThreshChi2Prob->push_back(chi2Prob);
		vector_twrs_emThreshTime->push_back(emTime);
		vector_twrs_emThreshRecoFlag->push_back(recoFlag);
		vector_twrs_emThreshSevLvl->push_back(sevlvl);
		vector_twrs_emThresh->push_back(emThresh);
		vector_twrs_emThreshEta->push_back(emThreshEta);
		vector_twrs_em3x3->push_back(em3x3);
		vector_twrs_em5x5->push_back(em5x5);
		vector_twrs_emSwiss->push_back(emSwiss);

		vector_twrs_hcalHitFlag->push_back( hcalFlag );
		vector_twrs_hcalHitDepth->push_back( hcalDepth );
		vector_twrs_hcalHitTime->push_back( hcalTime );
	}

	// put results into the event
	std::string branchprefix = aliasprefix_;
	if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

	iEvent.put(evt_ntwrs, "evtn"+branchprefix);
	iEvent.put(vector_twrs_eta, branchprefix+"eta");
	iEvent.put(vector_twrs_phi, branchprefix+"phi");
	iEvent.put(vector_twrs_detid, branchprefix+"detid");
	iEvent.put(vector_twrs_emEnergy, branchprefix+"emEnergy");
	iEvent.put(vector_twrs_hadEnergy, branchprefix+"hadEnergy");
	iEvent.put(vector_twrs_outerEnergy, branchprefix+"outerEnergy");

	iEvent.put(vector_twrs_emEt, branchprefix+"emEt");
	iEvent.put(vector_twrs_hadEt, branchprefix+"hadEt");
	iEvent.put(vector_twrs_outerEt, branchprefix+"outerEt");

	iEvent.put(vector_twrs_pcorr, branchprefix+"pcorr");
	iEvent.put(vector_twrs_etcorr, branchprefix+"etcorr");
	iEvent.put(vector_twrs_emEtcorr, branchprefix+"emEtcorr");
	iEvent.put(vector_twrs_hadEtcorr, branchprefix+"hadEtcorr");
	iEvent.put(vector_twrs_outerEtcorr, branchprefix+"outerEtcorr");
	iEvent.put(vector_twrs_etacorr, branchprefix+"etacorr");
	iEvent.put(vector_twrs_phicorr, branchprefix+"phicorr");

	iEvent.put(vector_twrs_ecalTime, branchprefix+"ecalTime");
	iEvent.put(vector_twrs_hcalTime, branchprefix+"hcalTime");
	iEvent.put(vector_twrs_hcalHitFlag, branchprefix+"hcalHitFlag");
	iEvent.put(vector_twrs_hcalHitTime, branchprefix+"hcalHitTime");
	iEvent.put(vector_twrs_hcalHitDepth, branchprefix+"hcalHitDepth");

	iEvent.put(vector_twrs_numBadEcalCells, branchprefix+"numBadEcalCells");
	iEvent.put(vector_twrs_numRecoveredEcalCells, branchprefix+"numRecoveredEcalCells");
	iEvent.put(vector_twrs_numProblematicEcalCells, branchprefix+"numProblematicEcalCells");

	iEvent.put(vector_twrs_numBadHcalCells, branchprefix+"numBadHcalCells");
	iEvent.put(vector_twrs_numRecoveredHcalCells, branchprefix+"numRecoveredHcalCells");
	iEvent.put(vector_twrs_numProblematicHcalCells, branchprefix+"numProblematicHcalCells");

	iEvent.put(vector_twrs_emMaxEcalMGPASampleADC, branchprefix+"emMaxEcalMGPASampleADC");
	iEvent.put(vector_twrs_numCrystals, branchprefix+"numCrystals");

	iEvent.put(vector_twrs_emThreshChi2, branchprefix+"emThreshChi2");
	//iEvent.put(vector_twrs_emThreshChi2Prob, branchprefix+"emThreshChi2Prob");
	iEvent.put(vector_twrs_emThreshTime, branchprefix+"emThreshTime");
	iEvent.put(vector_twrs_emThreshSevLvl, branchprefix+"emThreshSevLvl");
	iEvent.put(vector_twrs_emThreshRecoFlag, branchprefix+"emThreshRecoFlag");
	iEvent.put(vector_twrs_emThresh, branchprefix+"emThresh");
	iEvent.put(vector_twrs_emThreshEta, branchprefix+"emThreshEta");
	iEvent.put(vector_twrs_em3x3, branchprefix+"em3x3");
	iEvent.put(vector_twrs_em5x5, branchprefix+"em5x5");
	iEvent.put(vector_twrs_emSwiss, branchprefix+"emSwiss");

	iEvent.put(vector_twrs_spikeEt, branchprefix+"spikeEt");
	iEvent.put(vector_twrs_spikeR4, branchprefix+"spikeR4");

}


float CaloTowerMaker::recHitChi2(DetId emMaxId, const EcalRecHitCollection *recHits)
{               
        EcalRecHitCollection::const_iterator it = recHits->find(emMaxId);
        if (it != recHits->end()) return it->chi2();
        return -9999.99;
}

//in 3_5_5, the chi2Prob method will crash--cannot use. Replace with above.
float CaloTowerMaker::recHitChi2Prob(DetId emMaxId, const EcalRecHitCollection *recHits)
{               
        EcalRecHitCollection::const_iterator it = recHits->find(emMaxId);
        if (it != recHits->end()) return it->chi2Prob();
        return -9999.99;
}

float CaloTowerMaker::recHitTime(DetId emMaxId, const EcalRecHitCollection *recHits)
{
	EcalRecHitCollection::const_iterator it = recHits->find(emMaxId);
	if (it != recHits->end()) return it->time();
	return -9999.99;
}

int CaloTowerMaker::recHitFlag(DetId emMaxId, const EcalRecHitCollection *recHits)
{
	EcalRecHitCollection::const_iterator it = recHits->find(emMaxId);
	if (it != recHits->end()) return it->recoFlag();
	return -1;
}

int CaloTowerMaker::recHitSeverityLevel(DetId emMaxId, const EcalRecHitCollection *recHits)
{
	EcalRecHitCollection::const_iterator it = recHits->find(emMaxId);
	if (it != recHits->end()) {
	  const EcalSeverityLevelAlgo* theEcalSevLvlAlgo;
	  //const EcalChannelStatus* theEcalChStatus = new EcalChannelStatus(); //have to set this up, idiot
	  return theEcalSevLvlAlgo->severityLevel( emMaxId, *recHits, *theEcalChStatus_);
	}
	return -1;
}

void CaloTowerMaker::recHitSamples(DetId emMaxId, const EcalDigiCollection *digis, std::vector<int> &samples)
{


	samples.clear();
	if (digi_) {
		EcalDigiCollection::const_iterator it = digis->find(emMaxId);
		EcalDataFrame frame;
		if (it != digis->end()) {
			frame = (*it);
			for (size_t i = 0; i < 10; ++i) samples.push_back(frame.sample(i).adc());
		}
	}
	else {
		for (size_t i = 0; i < 10; ++i) samples.push_back(0);
	}

}


//returns sum of energies in emMax + 4 immediate neighbors
float CaloTowerMaker::SwissCross( reco::BasicCluster& dummyCluster, const EcalRecHitCollection *&recHits, const DetId& emMaxId) {
  EcalClusterTools clusterTools;
  float emSwiss = 0.;
  emSwiss += clusterTools.matrixEnergy(dummyCluster, recHits, topology_, emMaxId, 0, 0, -1, 1);  //topology_ is a data memeber 
  emSwiss += clusterTools.matrixEnergy(dummyCluster, recHits, topology_, emMaxId, -1, 1, 0, 0); 
  emSwiss -= clusterTools.matrixEnergy(dummyCluster, recHits, topology_, emMaxId, 0, 0, 0, 0); //center of cross was included twice above 
  return emSwiss;
}


// ------------ method called once each job just before starting event loop  ------------
	void 
CaloTowerMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
CaloTowerMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloTowerMaker);

