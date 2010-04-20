// -*- C++ -*-
//
// Package:    SCMaker
// Class:      SCMaker
// 
/**\class SCMaker SCMaker.cc CMS2/SCMaker/src/SCMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
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

#include "CMS2/NtupleMaker/interface/SCMaker.h"

#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

//
// class decleration
//

//
// constructors and destructor
//
SCMaker::SCMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // number of superclusters in the event
  produces<unsigned int>("evtnscs").setBranchAlias("evt_nscs");
  produces<unsigned int>("evtecalrecostatus").setBranchAlias("evt_ecalRecoStatus");

  // number of basicclusters and crystals
  produces<std::vector<float> >(branchprefix+"clustersSize").setBranchAlias(aliasprefix_+"_clustersSize");
  produces<std::vector<float> >(branchprefix+"crystalsSize").setBranchAlias(aliasprefix_+"_crystalsSize");

  // energies
  produces<std::vector<float> >(branchprefix+"energy").setBranchAlias(aliasprefix_+"_energy");
  produces<std::vector<float> >(branchprefix+"rawEnergy").setBranchAlias(aliasprefix_+"_rawEnergy"); 
  produces<std::vector<float> >(branchprefix+"preshowerEnergy").setBranchAlias(aliasprefix_+"_preshowerEnergy");

  // positions
  produces<std::vector<LorentzVector> >(branchprefix+"p4").setBranchAlias(aliasprefix_+"_p4");
  produces<std::vector<LorentzVector> >(branchprefix+"vtxp4").setBranchAlias(aliasprefix_+"_vtx_p4");
  produces<std::vector<LorentzVector> >(branchprefix+"posp4").setBranchAlias(aliasprefix_+"_pos_p4");
  produces<std::vector<float> >(branchprefix+"eta").setBranchAlias(aliasprefix_+"_eta");
  produces<std::vector<float> >(branchprefix+"phi").setBranchAlias(aliasprefix_+"_phi");

  // longitudinal shower shape and hcal isolations
  produces<std::vector<float> >(branchprefix+"hoe").setBranchAlias(aliasprefix_+"_hoe");
  //produces<std::vector<float> >(branchprefix+"hd1").setBranchAlias(aliasprefix_+"_hd1");
  //produces<std::vector<float> >(branchprefix+"hd2").setBranchAlias(aliasprefix_+"_hd2");

  // shape variables for seed basiccluster
  // see
  // RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
  // revision=1.7

  produces<std::vector<int> >(branchprefix+"detIdSeed").setBranchAlias(aliasprefix_+"_detIdSeed");
  produces<std::vector<int> >(branchprefix+"severitySeed").setBranchAlias(aliasprefix_+"_severitySeed");
  produces<std::vector<float> >(branchprefix+"timeSeed").setBranchAlias(aliasprefix_+"_timeSeed");
  produces<std::vector<float> >(branchprefix+"eSeed").setBranchAlias(aliasprefix_+"_eSeed");
  produces<std::vector<float> >(branchprefix+"eMax").setBranchAlias(aliasprefix_+"_eMax");
  produces<std::vector<float> >(branchprefix+"e2nd").setBranchAlias(aliasprefix_+"_e2nd");

  produces<std::vector<float> >(branchprefix+"e1x3").setBranchAlias(aliasprefix_+"_e1x3");
  produces<std::vector<float> >(branchprefix+"e3x1").setBranchAlias(aliasprefix_+"_e3x1"); 
  produces<std::vector<float> >(branchprefix+"e1x5").setBranchAlias(aliasprefix_+"_e1x5");
  produces<std::vector<float> >(branchprefix+"e2x2").setBranchAlias(aliasprefix_+"_e2x2"); 
  produces<std::vector<float> >(branchprefix+"e3x2").setBranchAlias(aliasprefix_+"_e3x2"); 
  produces<std::vector<float> >(branchprefix+"e3x3").setBranchAlias(aliasprefix_+"_e3x3"); 
  produces<std::vector<float> >(branchprefix+"e4x4").setBranchAlias(aliasprefix_+"_e4x4"); 
  produces<std::vector<float> >(branchprefix+"e5x5").setBranchAlias(aliasprefix_+"_e5x5"); 
  produces<std::vector<float> >(branchprefix+"e2x5Max").setBranchAlias(aliasprefix_+"_e2x5Max");
  // covariances
  produces<std::vector<float> >(branchprefix+"sigmaEtaEta").setBranchAlias(aliasprefix_+"_sigmaEtaEta");
  produces<std::vector<float> >(branchprefix+"sigmaEtaPhi").setBranchAlias(aliasprefix_+"_sigmaEtaPhi");
  produces<std::vector<float> >(branchprefix+"sigmaPhiPhi").setBranchAlias(aliasprefix_+"_sigmaPhiPhi");
  produces<std::vector<float> >(branchprefix+"sigmaIEtaIEta").setBranchAlias(aliasprefix_+"_sigmaIEtaIEta");
  produces<std::vector<float> >(branchprefix+"sigmaIEtaIPhi").setBranchAlias(aliasprefix_+"_sigmaIEtaIPhi");
  produces<std::vector<float> >(branchprefix+"sigmaIPhiIPhi").setBranchAlias(aliasprefix_+"_sigmaIPhiIPhi");

  // match to electrons
  produces<std::vector<int> >(branchprefix+"elsidx").setBranchAlias(aliasprefix_+"_elsidx");

  debug_ = iConfig.getParameter<bool>("debug");
  if (debug_) {
    produces<std::vector<float>         >(branchprefix+"mcdr"            ).setBranchAlias(aliasprefix_+"_mc_dr"           );
    produces<std::vector<float>         >(branchprefix+"mcenergy"           ).setBranchAlias(aliasprefix_+"_mc_energy"          );
  }

  // add superclusters to the ntuple if they have ET > scEtMin_
  scEtMin_ = iConfig.getParameter<double>("scEtMin");

  // hcal depth isolation
  //isoExtRadius_ = iConfig.getParameter<double> ("isoExtRadius");
  //isoIntRadius_ = iConfig.getParameter<double> ("isoIntRadius");
  //isoEtMin_ = iConfig.getParameter<double> ("isoEtMin");

  // input tags for superclusters
  scInputTag_EE_ = iConfig.getParameter<edm::InputTag>("scInputTag_EE");
  scInputTag_EB_ = iConfig.getParameter<edm::InputTag>("scInputTag_EB");
  scInputTags_.clear();
  scInputTags_.push_back(scInputTag_EE_);
  scInputTags_.push_back(scInputTag_EB_);

  hitInputTags_.clear();
  ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
  ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");
  hitInputTags_.push_back(ecalRecHitsInputTag_EE_);
  hitInputTags_.push_back(ecalRecHitsInputTag_EB_);

  // other input tags
  hcalRecHitsInputTag_HBHE_ = iConfig.getParameter<edm::InputTag>("hcalRecHitsInputTag_HBHE");
  primaryVertexInputTag_ = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
  electronsInputTag_ = iConfig.getParameter<edm::InputTag>("electronsInputTag");
  //caloTowersInputTag_ = iConfig.getParameter<edm::InputTag>("caloTowersInputTag");

  MCTruthCollection_ = iConfig.getParameter<edm::InputTag>("MCTruthCollection");

  // initialise this
  cachedCaloGeometryID_ = 0;

}

void SCMaker::beginRun(edm::Run& iRun, const edm::EventSetup& iSetup)
{
  /*
    edm::ESHandle<EcalChannelStatus> chStatus;
    iSetup.get<EcalChannelStatusRcd>().get(chStatus);
    // where const EcalChannelStatusCode * channelStatus;
    channelStatus_ = chStatus.product();
  */
}

void SCMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get the calo geometry
  if (cachedCaloGeometryID_ != iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
    cachedCaloGeometryID_ = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
    iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
  }

  // get hcal rechits
  edm::Handle<HBHERecHitCollection> hcalRecHitsHandle;
  try {
    iEvent.getByLabel(hcalRecHitsInputTag_HBHE_, hcalRecHitsHandle);
  }
  catch ( cms::Exception& ex ) {
    edm::LogError("SCMakerError") << "Error! can't get the HCAL Hits";
  }

  // get hcal rechit metacollection 
  HBHERecHitMetaCollection *mhbhe = 0;
  if (!hcalRecHitsHandle.failedToGet()) {
    mhbhe =  new HBHERecHitMetaCollection(*hcalRecHitsHandle);
  }

  // get the primary vertices
  edm::Handle<reco::VertexCollection> vertexHandle;
  try {
    iEvent.getByLabel(primaryVertexInputTag_, vertexHandle);
  }
  catch ( cms::Exception& ex ) {
    edm::LogError("SCMakerError") << "Error! can't get the primary vertex";
  }
  const reco::VertexCollection *vertexCollection = vertexHandle.product();
  Point pv(0.0, 0.0, 0.0);
  if (vertexCollection->size() > 0) {
    pv = vertexCollection->at(0).position();
  }

  // get the electrons (for matching)
  edm::Handle<reco::GsfElectronCollection> electronsHandle;
  try {
    iEvent.getByLabel(electronsInputTag_, electronsHandle);
  }
  catch ( cms::Exception& ex ) {
    edm::LogError("SCMakerError") << "Error! can't get the electrons";
  }
  const reco::GsfElectronCollection *electronsCollection = electronsHandle.product();

  // get hoe variable
  HoECalculator hoeCalc(caloGeometry_);

  // ecal cluster shape variables
  EcalClusterLazyTools lazyTools(iEvent, iSetup,
				 ecalRecHitsInputTag_EB_, ecalRecHitsInputTag_EE_);

  // get hcal depth isolations
  //edm::Handle<CaloTowerCollection> caloTowersHandle;
  //iEvent.getByLabel(caloTowersInputTag_, caloTowersHandle);
  //const CaloTowerCollection *coloTowersCollection = caloTowersHandle.product();
  //EgammaTowerIsolation egammaIsoD1(isoExtRadius_, isoIntRadius_, isoEtMin_, 1, coloTowersCollection);
  //EgammaTowerIsolation egammaIsoD2(isoExtRadius_, isoIntRadius_, isoEtMin_, 2, coloTowersCollection);

  std::auto_ptr<unsigned int> evt_nscs (new unsigned int);
  std::auto_ptr<unsigned int> evt_ecalRecoStatus (new unsigned int);

  std::auto_ptr<std::vector<LorentzVector> > vector_scs_p4 (new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_scs_pos_p4 (new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_scs_vtx_p4 (new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<float> > vector_scs_eta (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_phi (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_clustersSize (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_crystalsSize (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_energy (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_preshowerEnergy (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_rawEnergy (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_hoe (new std::vector<float>);
  //std::auto_ptr<std::vector<float> > vector_scs_hd1 (new std::vector<float>);
  //std::auto_ptr<std::vector<float> > vector_scs_hd2 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_eMax (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_e2nd (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_eSeed (new std::vector<float>);
  std::auto_ptr<std::vector<int> > vector_scs_severitySeed (new std::vector<int>);
  std::auto_ptr<std::vector<int> > vector_scs_detIdSeed (new std::vector<int>);
  std::auto_ptr<std::vector<float> > vector_scs_timeSeed (new std::vector<float>);

  std::auto_ptr<std::vector<float> > vector_scs_e1x3 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_e3x1 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_e1x5 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_e2x2 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_e3x2 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_e3x3 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_e4x4 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_e5x5 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_e2x5Max (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_sigmaEtaEta (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_sigmaEtaPhi(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_sigmaPhiPhi(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_sigmaIEtaIEta (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_sigmaIEtaIPhi(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_scs_sigmaIPhiIPhi(new std::vector<float>);
  std::auto_ptr<std::vector<int> > vector_scs_elsidx(new std::vector<int>);

  // for debugging - will only be filled and put in event if debug_ == true
  std::auto_ptr<std::vector<float> > vector_scs_mc_energy (new std::vector<float>); 
  std::auto_ptr<std::vector<float> > vector_scs_mc_dr     (new std::vector<float>); 

  *evt_nscs = 0;
  *evt_ecalRecoStatus = 0;

  // there are multiple supercluster collections. In the ntuple
  // these will become concatonated
  for (unsigned int i = 0; i < scInputTags_.size(); ++i)
    {

      // get superclusters
      edm::Handle<reco::SuperClusterCollection> scHandle;
      try {
	iEvent.getByLabel(scInputTags_[i], scHandle);
      }
      catch ( cms::Exception& ex ) {
	edm::LogError("SCMakerError") << "Error! can't get the SuperClusters";
      }
      const reco::SuperClusterCollection *scCollection = scHandle.product();

      const HepMC::GenEvent* genEvent = 0;
      if (debug_) {
	edm::Handle<edm::HepMCProduct> pMCTruth ;
	iEvent.getByLabel(MCTruthCollection_, pMCTruth);
	genEvent = pMCTruth->GetEvent();
      }

      // get hits
      edm::Handle<EcalRecHitCollection> rhcHandle;
      iEvent.getByLabel(hitInputTags_[i], rhcHandle);
      const EcalRecHitCollection *recHits = rhcHandle.product();

      if (recHits->size() > 0) *evt_ecalRecoStatus |= (1<<i);

      size_t scIndex = 0;
      for (reco::SuperClusterCollection::const_iterator sc = scCollection->begin();
	   sc != scCollection->end(); ++sc, ++scIndex) {

	// do ET cut
	if ( (sc->energy()/cosh(sc->eta())) < scEtMin_) continue;

	if (debug_) {
	  // truth
	  double dRClosest = 999.9;
	  double energyClosest = 0;
	  closestMCParticle(genEvent, *sc, dRClosest, energyClosest);
	  // fill vector
	  vector_scs_mc_energy      ->push_back(energyClosest        );
	  vector_scs_mc_dr       ->push_back( dRClosest         );
	}

	LorentzVector p4 = initP4(pv, *sc);
	vector_scs_p4->push_back( p4 );
	vector_scs_vtx_p4->push_back( LorentzVector(pv.x(), pv.y(), pv.z(), 0.) );
	vector_scs_pos_p4->push_back( LorentzVector(sc->position().x(), sc->position().y(), sc->position().z(), 0.) );
	vector_scs_eta->push_back( sc->eta() );
	vector_scs_phi->push_back( sc->phi() );
	vector_scs_energy->push_back( sc->energy() );
	vector_scs_rawEnergy->push_back( sc->rawEnergy() );
	vector_scs_preshowerEnergy->push_back( sc->preshowerEnergy() );
	vector_scs_hoe->push_back( hoeCalc(&(*sc), mhbhe) );
	//vector_scs_hd1->push_back(egammaIsoD1.getTowerEtSum(&(*sc)) );
	//vector_scs_hd2->push_back(egammaIsoD2.getTowerEtSum(&(*sc)) );

	DetId seedId = sc->seed()->seed();
	EcalRecHitCollection::const_iterator seedHit = recHits->find(seedId);
	if (seedHit != recHits->end()) {
	  vector_scs_eSeed->push_back( seedHit->energy() );
	  vector_scs_detIdSeed->push_back(seedHit->id().rawId());
	  vector_scs_severitySeed->push_back ( seedHit->recoFlag() );
	  vector_scs_timeSeed->push_back (seedHit->time() );
	} else {
	  vector_scs_eSeed->push_back(0.0);
	  vector_scs_detIdSeed->push_back(-1);
	  vector_scs_severitySeed->push_back (-1);
	  vector_scs_timeSeed->push_back (-9999.99);
	}

	vector_scs_eMax->push_back( lazyTools.eMax(*(sc->seed())) );
	vector_scs_e2nd->push_back( lazyTools.e2nd(*(sc->seed())) );

	vector_scs_e1x3->push_back( lazyTools.e1x3(*(sc->seed())) );
	vector_scs_e3x1->push_back( lazyTools.e3x1(*(sc->seed())) );
	vector_scs_e1x5->push_back( lazyTools.e1x5(*(sc->seed())) );
	vector_scs_e2x2->push_back( lazyTools.e2x2(*(sc->seed())) );
	vector_scs_e3x2->push_back( lazyTools.e3x2(*(sc->seed())) );
	vector_scs_e3x3->push_back( lazyTools.e3x3(*(sc->seed())) );
	vector_scs_e4x4->push_back( lazyTools.e4x4(*(sc->seed())) );
	vector_scs_e5x5->push_back( lazyTools.e5x5(*(sc->seed())) );
	vector_scs_e2x5Max->push_back( lazyTools.e2x5Max(*(sc->seed())) );
	std::vector<float> covariances = lazyTools.covariances(*(sc->seed()));
	// if seed basic cluster is in the endcap then correct sigma eta eta
	// according to the super cluster eta
	if(fabs(sc->seed()->eta()) > 1.479) {
	  covariances[0] -= 0.02*(fabs(sc->eta()) - 2.3);
	}
	vector_scs_sigmaEtaEta->push_back( covariances[0] > 0 ? sqrt(covariances[0]) : -sqrt(-covariances[0]) );
	vector_scs_sigmaEtaPhi->push_back( covariances[1] > 0 ? sqrt(covariances[1]) : -sqrt(-covariances[1]) );
	vector_scs_sigmaPhiPhi->push_back( covariances[2] > 0 ? sqrt(covariances[2]) : -sqrt(-covariances[2]) );
	std::vector<float> localCovariances = lazyTools.localCovariances(*(sc->seed()));
	vector_scs_sigmaIEtaIEta->push_back( localCovariances[0] > 0 ? sqrt(localCovariances[0]) : -sqrt(-localCovariances[0]) );
	vector_scs_sigmaIEtaIPhi->push_back( localCovariances[1] > 0 ? sqrt(localCovariances[1]) : -sqrt(-localCovariances[1]) );
	vector_scs_sigmaIPhiIPhi->push_back( localCovariances[2] > 0 ? sqrt(localCovariances[2]) : -sqrt(-localCovariances[2]) );
	vector_scs_clustersSize->push_back( sc->clustersSize() );
	const std::vector<std::pair<DetId, float > > detIds = sc->hitsAndFractions() ;
	vector_scs_crystalsSize->push_back( detIds.size() );

	// do match to electrons
	const edm::Ref<reco::SuperClusterCollection> scRef(scHandle, scIndex);
	int electronIndex = -9999;
	for (size_t i = 0; i < electronsCollection->size(); ++i)
	  {
	    if ((*electronsCollection)[i].superCluster() == scRef) {
	      electronIndex = i;
	      break;
	    }
	  }	
	vector_scs_elsidx->push_back(electronIndex);

      } // end loop on scs

    } // end loop on sc input tags

  *evt_nscs = vector_scs_p4->size();

  // put results into the event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(evt_nscs, "evtnscs");
  iEvent.put(evt_ecalRecoStatus, "evtecalrecostatus");

  iEvent.put(vector_scs_energy, branchprefix+"energy");
  iEvent.put(vector_scs_rawEnergy, branchprefix+"rawEnergy");
  iEvent.put(vector_scs_preshowerEnergy, branchprefix+"preshowerEnergy");
  iEvent.put(vector_scs_p4, branchprefix+"p4");
  iEvent.put(vector_scs_vtx_p4, branchprefix+"vtxp4");
  iEvent.put(vector_scs_pos_p4, branchprefix+"posp4");
  iEvent.put(vector_scs_eta, branchprefix+"eta");
  iEvent.put(vector_scs_phi, branchprefix+"phi");
  iEvent.put(vector_scs_hoe, branchprefix+"hoe");
  //iEvent.put(vector_scs_hd1, branchprefix+"hd1");
  //iEvent.put(vector_scs_hd2, branchprefix+"hd2");
  iEvent.put(vector_scs_eSeed, branchprefix+"eSeed");
  iEvent.put(vector_scs_detIdSeed, branchprefix+"detIdSeed");
  iEvent.put(vector_scs_severitySeed, branchprefix+"severitySeed");
  iEvent.put(vector_scs_timeSeed, branchprefix+"timeSeed");
  iEvent.put(vector_scs_e2nd, branchprefix+"e2nd");
  iEvent.put(vector_scs_eMax, branchprefix+"eMax");

  iEvent.put(vector_scs_e1x3, branchprefix+"e1x3");
  iEvent.put(vector_scs_e3x1, branchprefix+"e3x1");
  iEvent.put(vector_scs_e1x5, branchprefix+"e1x5");
  iEvent.put(vector_scs_e2x2, branchprefix+"e2x2");
  iEvent.put(vector_scs_e3x2, branchprefix+"e3x2");
  iEvent.put(vector_scs_e3x3, branchprefix+"e3x3");
  iEvent.put(vector_scs_e4x4, branchprefix+"e4x4");
  iEvent.put(vector_scs_e5x5, branchprefix+"e5x5");
  iEvent.put(vector_scs_e2x5Max, branchprefix+"e2x5Max");
  iEvent.put(vector_scs_sigmaEtaEta, branchprefix+"sigmaEtaEta");
  iEvent.put(vector_scs_sigmaEtaPhi, branchprefix+"sigmaEtaPhi");
  iEvent.put(vector_scs_sigmaPhiPhi, branchprefix+"sigmaPhiPhi");
  iEvent.put(vector_scs_sigmaIEtaIEta, branchprefix+"sigmaIEtaIEta");
  iEvent.put(vector_scs_sigmaIEtaIPhi, branchprefix+"sigmaIEtaIPhi");
  iEvent.put(vector_scs_sigmaIPhiIPhi, branchprefix+"sigmaIPhiIPhi");
  iEvent.put(vector_scs_clustersSize, branchprefix+"clustersSize");
  iEvent.put(vector_scs_crystalsSize, branchprefix+"crystalsSize");
  iEvent.put(vector_scs_elsidx, branchprefix+"elsidx");

  if (debug_) {
    iEvent.put(vector_scs_mc_dr           ,branchprefix+"mcdr"          );
    iEvent.put(vector_scs_mc_energy          ,branchprefix+"mcenergy"         );
  }

  delete mhbhe;

}

math::XYZTLorentzVectorF SCMaker::initP4(const math::XYZPoint &pvPos, 
					 const reco::SuperCluster &sc)
{

  math::XYZVector scPos(sc.x(), sc.y(), sc.z());
  math::XYZVector pvPosVec(pvPos.x(), pvPos.y(), pvPos.z());
  math::XYZVector objPosition = scPos - pvPosVec;
  double scale = sc.energy() / objPosition.R();
  return math::XYZTLorentzVectorF(objPosition.x() * scale, 
				  objPosition.y() * scale, 
				  objPosition.z() * scale, 
				  sc.energy());
}



// ------------ method called once each job just before starting event loop  ------------
void 
SCMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCMaker::endJob() 
{
}

//
// Closest MC Particle
//
void SCMaker::closestMCParticle(const HepMC::GenEvent *genEvent, const reco::SuperCluster &sc,
				double &dRClosest, double &energyClosest)
{

  // SuperCluster eta, phi
  double scEta = sc.eta();
  double scPhi = sc.phi();

  // initialize dRClosest to a large number
  dRClosest = 999.9;

  // loop over the MC truth particles to find the
  // closest to the superCluster in dR space
  for(HepMC::GenEvent::particle_const_iterator currentParticle = genEvent->particles_begin();
      currentParticle != genEvent->particles_end(); currentParticle++ )
    {
      if((*currentParticle)->status() == 3 && abs((*currentParticle)->pdg_id()) == 11)
	{
	  // need GenParticle in ECAL co-ordinates
	  HepMC::FourVector vtx = (*currentParticle)->production_vertex()->position();
	  double phiTrue = (*currentParticle)->momentum().phi();
	  double etaTrue = ecalEta((*currentParticle)->momentum().eta(), vtx.z()/10., vtx.perp()/10.);

	  double dPhi = reco::deltaPhi(phiTrue, scPhi);
	  double dEta = scEta - etaTrue;
	  double deltaR = std::sqrt(dPhi*dPhi + dEta*dEta);

	  if(deltaR < dRClosest)
	    {
	      dRClosest = deltaR;
	      energyClosest = (*currentParticle)->momentum().e();
	    }

	} // end if stable particle     

    } // end loop on get particles

}

//
// Compute Eta in the ECAL co-ordinate system
//
float SCMaker::ecalEta(float EtaParticle , float Zvertex, float plane_Radius)
{
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479;

  if(EtaParticle != 0.)
    {
      float Theta = 0.0  ;
      float ZEcal = (R_ECAL-plane_Radius)*sinh(EtaParticle)+Zvertex;

      if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
      if(Theta<0.0) Theta = Theta+Geom::pi() ;

      float ETA = - log(tan(0.5*Theta));

      if( fabs(ETA) > etaBarrelEndcap )
	{
	  float Zend = Z_Endcap ;
	  if(EtaParticle<0.0 )  Zend = -Zend ;
	  float Zlen = Zend - Zvertex ;
	  float RR = Zlen/sinh(EtaParticle);
	  Theta = atan((RR+plane_Radius)/Zend);
	  if(Theta<0.0) Theta = Theta+Geom::pi() ;
	  ETA = - log(tan(0.5*Theta));
	}

      return ETA;
    }
  else
    {
      edm::LogWarning("")  << "[EgammaSuperClusters::ecalEta] Warning: Eta equals to zero, not correcting" ;
      return EtaParticle;
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SCMaker);

