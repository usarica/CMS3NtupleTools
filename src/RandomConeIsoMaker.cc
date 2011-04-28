//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      RandomConeIsoMaker
// 
/**\class RandomConeIsoMaker RandomConeIsoMaker.cc CMS2/NtupleMakerMaker/src/RandomConeIsoMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: RandomConeIsoMaker.cc,v 1.7 2011/04/28 00:59:56 dbarge Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"

// random numver generation
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FastSimulation/Utilities/interface/RandomEngine.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CLHEP/Random/JamesRandom.h"

// propagation to calorimeter
#include "PhysicsTools/IsolationAlgos/interface/PropagateToCal.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfo.h"

// geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Vector3D.h"


// reco objects
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EcalScDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDigi/interface/EcalSrFlag.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDigi/interface/EBSrFlag.h"
#include "DataFormats/EcalDigi/interface/EESrFlag.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// root
#include "Math/VectorUtil.h"
#include "TVector3.h"
// isolations
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractorFactory.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"

#include "CMS2/NtupleMaker/interface/RandomConeIsoMaker.h"




typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

RandomConeIsoMaker::RandomConeIsoMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");


  produces<std::vector<float> >               (branchprefix+"ecalIso03egamma"         ).setBranchAlias(aliasprefix_+"_ecalIso03_egamma"        );
  produces<std::vector<float> >               (branchprefix+"hcalIso03egamma"         ).setBranchAlias(aliasprefix_+"_hcalIso03_egamma"        );
  produces<std::vector<float> >               (branchprefix+"hcalD1Iso03egamma"       ).setBranchAlias(aliasprefix_+"_hcalD1Iso03_egamma"      );
  produces<std::vector<float> >               (branchprefix+"hcalD2Iso03egamma"       ).setBranchAlias(aliasprefix_+"_hcalD2Iso03_egamma"      );
  produces<std::vector<float> >               (branchprefix+"trkIso03egamma"          ).setBranchAlias(aliasprefix_+"_trkIso03_egamma"         );
  
  produces<std::vector<float> >               (branchprefix+"ecalIso03mu"             ).setBranchAlias(aliasprefix_+"_ecalIso03_mu"            );
  produces<std::vector<float> >               (branchprefix+"hcalIso03mu"             ).setBranchAlias(aliasprefix_+"_hcalIso03_mu"            );
  produces<std::vector<float> >               (branchprefix+"hoIso03mu"               ).setBranchAlias(aliasprefix_+"_hoIso03_mu"              );
  produces<std::vector<float> >               (branchprefix+"trkIso03mu"              ).setBranchAlias(aliasprefix_+"_trkIso03_mu"             );
 
  produces<std::vector<int> >                 (branchprefix+"srflag"                  ).setBranchAlias(aliasprefix_+"_srFlag"                  );
  produces<std::vector<float> >               (branchprefix+"toweremet"               ).setBranchAlias(aliasprefix_+"_towerEmEt"               );
  produces<std::vector<float> >               (branchprefix+"towerhadet"              ).setBranchAlias(aliasprefix_+"_towerHadEt"              );
  produces<std::vector<float> >               (branchprefix+"drclosesttower"          ).setBranchAlias(aliasprefix_+"_dRClosestTower"          );
  produces<std::vector<float> >               (branchprefix+"drclosesttoweremet"      ).setBranchAlias(aliasprefix_+"_dRClosestTowerEmEt"      );
  
  produces<std::vector<LorentzVector> >       (branchprefix+"trksp4"                  ).setBranchAlias(aliasprefix_+"_trksp4"                  );
  produces<std::vector<LorentzVector> >       (branchprefix+"trksecalp4"              ).setBranchAlias(aliasprefix_+"_trksecalp4"              );

  
  //
  // Set up random number distribution
  // See Example 1 of 
  // SWGuideEDMRandomNumberGeneratorService
  // 
  
  edm::Service<edm::RandomNumberGenerator> rng;
  if (!rng.isAvailable())
    throw cms::Exception("Configuration") << "RandomNumberGeneratorService not present\n";
  std::cout << "[EgammaIsolationAnalyzer::EgammaIsolationAnalyzer] seed: " << rng->mySeed() << std::endl;
  jamesRandom_ = new CLHEP::HepJamesRandom(rng->mySeed());
  
  primaryVertexInputTag_ = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
  
  
  ecalBarrelRecHitProducer_       = iConfig.getParameter<edm::InputTag>("ecalBarrelRecHitProducer");
  ecalBarrelRecHitCollection_     = iConfig.getParameter<edm::InputTag>("ecalBarrelRecHitCollection");
  ecalEndcapRecHitProducer_       = iConfig.getParameter<edm::InputTag>("ecalEndcapRecHitProducer");
  ecalEndcapRecHitCollection_     = iConfig.getParameter<edm::InputTag>("ecalEndcapRecHitCollection");
  trackProducer_                  = iConfig.getParameter<edm::InputTag>("trackProducer");
  beamspotProducer_               = iConfig.getParameter<edm::InputTag>("BeamspotProducer");
  towerProducer_                  = iConfig.getParameter<edm::InputTag>("towerProducer");


   // Load MuIsoExtractor parameters
  edm::ParameterSet caloExtractorPSet = iConfig.getParameter<edm::ParameterSet>("CaloExtractorPSet");
  std::string caloExtractorName = caloExtractorPSet.getParameter<std::string>("ComponentName");
  muIsoExtractorCalo_ = IsoDepositExtractorFactory::get()->create( caloExtractorName, caloExtractorPSet);
 
  edm::ParameterSet trackExtractorPSet = iConfig.getParameter<edm::ParameterSet>("TrackExtractorPSet");
  std::string trackExtractorName = trackExtractorPSet.getParameter<std::string>("ComponentName");
  muIsoExtractorTrack_ = IsoDepositExtractorFactory::get()->create( trackExtractorName, trackExtractorPSet);
  
  edm::ParameterSet jetExtractorPSet = iConfig.getParameter<edm::ParameterSet>("JetExtractorPSet");
  std::string jetExtractorName = jetExtractorPSet.getParameter<std::string>("ComponentName");
  muIsoExtractorJet_ = IsoDepositExtractorFactory::get()->create( jetExtractorName, jetExtractorPSet);

// Load egamma Iso parameters
  edm::ParameterSet    	egammaPSet = iConfig.getParameter<edm::ParameterSet>("egammaPSet");
  
  //vetos
  egIsoPtMinBarrel_               = egammaPSet.getParameter<double>("etMinBarrel");
  egIsoEMinBarrel_                = egammaPSet.getParameter<double>("eMinBarrel");
  egIsoPtMinEndcap_               = egammaPSet.getParameter<double>("etMinEndcap");
  egIsoEMinEndcap_                = egammaPSet.getParameter<double>("eMinEndcap");
  egIsoConeSizeInBarrel_          = egammaPSet.getParameter<double>("intRadiusBarrel");
  egIsoConeSizeInEndcap_          = egammaPSet.getParameter<double>("intRadiusEndcap");
  egIsoConeSizeOut_               = egammaPSet.getParameter<double>("extRadius");
  egIsoJurassicWidth_             = egammaPSet.getParameter<double>("jurassicWidth");
  
  
  
  // options
  useIsolEt_      = egammaPSet.getParameter<bool>("useIsolEt");
  tryBoth_        = egammaPSet.getParameter<bool>("tryBoth");
  subtract_       = egammaPSet.getParameter<bool>("subtract");
  useNumCrystals_ = egammaPSet.getParameter<bool>("useNumCrystals");
  vetoClustered_  = egammaPSet.getParameter<bool>("vetoClustered");
  
  
  ////hcal in egamma

 
  egHcalIsoPtMin_               = egammaPSet.getParameter<double>("etMin");
  egHcalIsoConeSizeIn_          = egammaPSet.getParameter<double>("intRadius");
  egHcalIsoConeSizeOut_         = egammaPSet.getParameter<double>("extRadius");
  egHcalDepth_                  = egammaPSet.getParameter<int>("Depth");
  
  
  ///trkiso in egamma

  ptMin_                = egammaPSet.getParameter<double>("ptMin_trk");
  intRadius_            = egammaPSet.getParameter<double>("intRadius_trk");
  extRadius_            = egammaPSet.getParameter<double>("extRadius_trk");
  maxVtxDist_           = egammaPSet.getParameter<double>("maxVtxDist_trk");
  drb_                  = egammaPSet.getParameter<double>("maxVtxDistXY_trk");
  
  	// SR
	srProducerEE_		= iConfig.getParameter<edm::InputTag>("srProducerEE");
	srProducerEB_           = iConfig.getParameter<edm::InputTag>("srProducerEB");

	// TrackAssociator parameters
	edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
	parameters_.loadParameters( parameters );
	trackAssociator_.useDefaultPropagator();


  // datasetName_ = iConfig.getParameter<std::string>("datasetName");
//   CMS2tag_     = iConfig.getParameter<std::string>("CMS2tag");
}


RandomConeIsoMaker::~RandomConeIsoMaker() {}

void RandomConeIsoMaker::beginJob() {  
}

void RandomConeIsoMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void RandomConeIsoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  

  std::auto_ptr<std::vector<float> >    ran_ecalIso03_egamma       (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_hcalIso03_egamma       (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_hcalD1Iso03_egamma     (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_hcalD2Iso03_egamma     (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_trkIso03_egamma        (new  std::vector<float>              );
  
  std::auto_ptr<std::vector<float> >    ran_ecalIso03_mu           (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_hcalIso03_mu           (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_hoIso03_mu             (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_trkIso03_mu            (new  std::vector<float>              );

  std::auto_ptr<std::vector<int> >      ran_srFlag                 (new  std::vector<int>                );
  std::auto_ptr<std::vector<float> >    ran_towerEmEt              (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_towerHadEt             (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_dRClosestTower         (new  std::vector<float>              );
  std::auto_ptr<std::vector<float> >    ran_dRClosestTowerEmEt     (new  std::vector<float>              );


  std::auto_ptr<std::vector<LorentzVector> > ran_trksp4            (new  std::vector<LorentzVector>      );
  std::auto_ptr<std::vector<LorentzVector> > ran_trksecalp4        (new  std::vector<LorentzVector>      );


    // get the primary vertices
  edm::Handle<reco::VertexCollection> vertexHandle;
  try {
    iEvent.getByLabel(primaryVertexInputTag_, vertexHandle);
  }
  catch ( cms::Exception& ex ) {
    edm::LogError("VertexMakerError") << "Error! can't get the primary vertex";
  }

  const reco::VertexCollection *vertexCollection = vertexHandle.product();
  //
  // get products
  //
  
  edm::Handle<EcalRecHitCollection> ecalBarrelRecHitHandle; //EcalRecHitCollection is a typedef to 
  iEvent.getByLabel(ecalBarrelRecHitProducer_.label(),ecalBarrelRecHitCollection_.label(), ecalBarrelRecHitHandle);
  EcalRecHitMetaCollection ecalBarrelHits(*ecalBarrelRecHitHandle);
  
  edm::Handle<EcalRecHitCollection> ecalEndcapRecHitHandle;
  iEvent.getByLabel(ecalEndcapRecHitProducer_.label(), ecalEndcapRecHitCollection_.label(),ecalEndcapRecHitHandle);
  EcalRecHitMetaCollection ecalEndcapHits(*ecalEndcapRecHitHandle);
  
  edm::Handle<CaloTowerCollection> towerHandle;
  iEvent.getByLabel(towerProducer_, towerHandle);
  const CaloTowerCollection* towers = towerHandle.product();
  std::map<DetId, CaloTower> towersMap;
  for (size_t i = 0; i < towers->size(); ++i) {
    towersMap.insert(std::make_pair<DetId, CaloTower>((*towers)[i].id(), (*towers)[i]));
  }
  
  edm::Handle<EESrFlagCollection> srHandleEE;
  iEvent.getByLabel(srProducerEE_, srHandleEE);
  const EESrFlagCollection *srFlagCollectionEE = 0;
  std::map<EcalScDetId, int> eeSrMap;
  
  if (iEvent.getByLabel(srProducerEE_, srHandleEE)) {
    srFlagCollectionEE = srHandleEE.product();
    for (size_t i = 0; i < srFlagCollectionEE->size(); ++i) {
      eeSrMap.insert(std::make_pair<EcalScDetId, int>((*srFlagCollectionEE)[i].id(), (*srFlagCollectionEE)[i].value()));
    }
  }
  
  edm::Handle<EBSrFlagCollection> srHandleEB;
  iEvent.getByLabel(srProducerEB_, srHandleEB);
  const EBSrFlagCollection *srFlagCollectionEB = 0;
  std::map<EcalTrigTowerDetId, int> ebSrMap;
  
  if (iEvent.getByLabel(srProducerEB_, srHandleEB)) {
    srFlagCollectionEB = srHandleEB.product();
    for (size_t i = 0; i < srFlagCollectionEB->size(); ++i) {
      ebSrMap.insert(std::make_pair<EcalTrigTowerDetId, int>((*srFlagCollectionEB)[i].id(), (*srFlagCollectionEB)[i].value()));
    }
  }
  
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackProducer_,tracks);
  const reco::TrackCollection* trackCollection = tracks.product();
  
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByLabel(beamspotProducer_,beamSpotH);
  reco::TrackBase::Point beamspot = beamSpotH->position();
  
  //Get Calo Geometry
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* caloGeom = pG.product();
  
  
  //
  // Set up the isolation algos
  //
  
  // Standard ecal isolation code
  
  EgammaRecHitIsolation ecalBarrelIsol( egIsoConeSizeOut_, egIsoConeSizeInBarrel_, egIsoJurassicWidth_, egIsoPtMinBarrel_, egIsoEMinBarrel_, caloGeom, &ecalBarrelHits, 0, DetId::Ecal );
  ecalBarrelIsol.setUseNumCrystals(useNumCrystals_);
  ecalBarrelIsol.setVetoClustered(vetoClustered_);
  
  EgammaRecHitIsolation ecalEndcapIsol( egIsoConeSizeOut_, egIsoConeSizeInEndcap_, egIsoJurassicWidth_, egIsoPtMinEndcap_, egIsoEMinEndcap_, caloGeom, &ecalEndcapHits, 0, DetId::Ecal);
  ecalEndcapIsol.setUseNumCrystals(useNumCrystals_);
  ecalEndcapIsol.setVetoClustered(vetoClustered_);
  
  // Standard hcal isolation code
  
  EgammaTowerIsolation myHadIsolationD1(egHcalIsoConeSizeOut_, egHcalIsoConeSizeIn_, egHcalIsoPtMin_, 1, towers);
  EgammaTowerIsolation myHadIsolationD2(egHcalIsoConeSizeOut_, egHcalIsoConeSizeIn_, egHcalIsoPtMin_, 2, towers);
  
  EgammaTowerIsolation hcalIsol(egHcalIsoConeSizeOut_, egHcalIsoConeSizeIn_, egHcalIsoPtMin_, egHcalDepth_, towers) ;
  //
  // define the properties of the direction of propagation
  // note: the mass and pt don't matter but is merely needed for convenience
  // in construction of objects to transform from one co-ordinate system
  // to another
  //
 
  GlobalPoint vertex(0.0, 0.0, 0.0);
  for (reco::VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx) {
    if(!vtx->isFake()){
      GlobalPoint vertex_tmp( vtx ->position().x(), vtx->position().y(),vtx->position().z());
      vertex = vertex_tmp;
    }
  }

   
  double mass = 0.000511;
  double pt = 20.0;
  int charge = 1;
  
 
  
  //
  // make a direction in eta-phi at random
  // and then propagate it to the calorimeter
  // such that an x-y-z position at the ecal is found
  // to construct a supercluster from
  //
  
  std::vector<std::pair<double, double> > randomDirections;
  double eta_eb = signedRnd(jamesRandom_->flat(), 0.0, 1.5);
  double phi = signedRnd(jamesRandom_->flat(), 0.0, 1.0) * M_PI;
  double eta_ee = signedRnd(jamesRandom_->flat(), 1.5, 3.0);
  randomDirections.push_back(std::make_pair<double, double>(eta_eb, phi));
  randomDirections.push_back(std::make_pair<double, double>(-1*eta_eb, phi));
  randomDirections.push_back(std::make_pair<double, double>(eta_ee, phi));
  randomDirections.push_back(std::make_pair<double, double>(-1*eta_ee, phi));

  for (size_t d = 0; d < randomDirections.size(); ++d)
    {
      
      math::PtEtaPhiMLorentzVector direction(pt, randomDirections[d].first, randomDirections[d].second, mass);
      math::XYZVector directionAtVertex(direction.px(), direction.py(), direction.pz());
      TrackDetMatchInfo info = trackAssociator_.associate(iEvent, iSetup, 
				GlobalVector(direction.px(), direction.py(), direction.pz()), 
							  vertex, charge, parameters_);
      
      // get the cone back to back with it
      // do once for EB and once for EE.
      
      
      const std::vector<DetId> &ecalIds = info.crossedEcalIds;
      const std::vector<const CaloTower *> &crossedTowers = info.crossedTowers;
      
      float hadEt = 0.0;
      float emEt = 0.0;
      if (crossedTowers.size() > 0) {
	hadEt = crossedTowers[0]->hadEt();
	emEt = crossedTowers[0]->emEt();
      }
      
      // find closest calo tower with Et > 1.0 GeV
      float dRClosest = 999;
      float etClosest = 0.0;
      for (size_t i = 0; i < towers->size(); ++i) {
	if ((*towers)[i].emEt() < 1.0) continue;
		  
	double dEta = (*towers)[i].eta() - info.trkGlobPosAtEcal.eta();
	double dPhi = acos(cos((*towers)[i].phi() - info.trkGlobPosAtEcal.phi()));
	double dR = sqrt(dEta*dEta + dPhi*dPhi);
	
	if (dR < dRClosest) {
	  dRClosest = dR;
	  etClosest = (*towers)[i].emEt();
	}
      }
      
      ran_dRClosestTowerEmEt->push_back(etClosest);
      ran_dRClosestTower->push_back(dRClosest);
      
      int srFlag = -1;
      if (ecalIds.size() > 0) {
	// EB
	if ((ecalIds[0].subdetId() == EcalBarrel) && srFlagCollectionEB != 0) {
	  EBDetId ebId(ecalIds[0]);
	  EcalTrigTowerDetId trigTowerId = ebId.tower();
	  std::map<EcalTrigTowerDetId, int>::const_iterator it = ebSrMap.find(trigTowerId);
	  if (it != ebSrMap.end()) srFlag = it->second;
	}
	// EE
	if ((ecalIds[0].subdetId() == EcalEndcap) && srFlagCollectionEE != 0) {
	  EEDetId eeId(ecalIds[0]);
	  const int scEdge = 5;
	  EcalScDetId trigTowerId((eeId.ix()-1)/scEdge+1, (eeId.iy()-1)/scEdge+1, eeId.zside());
	  std::map<EcalScDetId, int>::const_iterator it = eeSrMap.find(trigTowerId);
	  if (it != eeSrMap.end()) srFlag = it->second;
	}
      }
      
      ran_srFlag->push_back(srFlag);
      ran_towerEmEt->push_back(emEt);
      ran_towerHadEt->push_back(hadEt);
      
      // make a super cluster and get a ref to it
      const math::XYZPoint positionAtCalo(info.trkGlobPosAtEcal.x(), info.trkGlobPosAtEcal.y(), info.trkGlobPosAtEcal.z()); 
      reco::SuperCluster superCluster(pt, positionAtCalo);
      reco::SuperClusterCollection *temporarySuperClusterCollection = new reco::SuperClusterCollection;
      temporarySuperClusterCollection->push_back(superCluster);
      reco::SuperClusterRef tempSCRef(temporarySuperClusterCollection, 0);
      
      // make an ecal candidate
      reco::RecoEcalCandidate ecalCand;
      ecalCand.setSuperCluster(tempSCRef);
      
      //
      // now everything is ready to use the standard ecal 
      // isolation code
      //
      
      double isolation = 0.0;
      
      if(tryBoth_){ //barrel + endcap
			if(useIsolEt_) isolation =  ecalBarrelIsol.getEtSum(&ecalCand) + ecalEndcapIsol.getEtSum(&ecalCand);
			else           isolation =  ecalBarrelIsol.getEnergySum(&ecalCand) + ecalEndcapIsol.getEnergySum(&ecalCand);
      }
      else if(fabs(superCluster.eta()) < 1.479) { //barrel
	std::cout << "barrel" << std::endl;
	if(useIsolEt_) isolation =  ecalBarrelIsol.getEtSum(&ecalCand);
	else           isolation =  ecalBarrelIsol.getEnergySum(&ecalCand);
      }
      else{ //endcap
	std::cout << "endcap" << std::endl;
	if(useIsolEt_) isolation =  ecalEndcapIsol.getEtSum(&ecalCand);
	else           isolation =  ecalEndcapIsol.getEnergySum(&ecalCand);
      }
      
      ran_ecalIso03_egamma->push_back(isolation);
      
      //
      // Standard hcal isolation code
      //

      double isolationD1 = myHadIsolationD1.getTowerEtSum(&superCluster);
      double isolationD2 = myHadIsolationD2.getTowerEtSum(&superCluster);
      double isolationhcal = hcalIsol.getTowerEtSum(&superCluster);
      
      ran_hcalD1Iso03_egamma->push_back(isolationD1);
      ran_hcalD2Iso03_egamma->push_back(isolationD2);
      ran_hcalIso03_egamma  ->push_back(isolationhcal);
      
      //
      // standard track isolation code
      //
      
      std::pair<int, double> isolationTk = getTrackIso(beamspot, trackCollection, directionAtVertex);  
      ran_trkIso03_egamma->push_back(isolationTk.second);

      ////standard muon isolation code
      TVector3 v3(0,0,0);
      double chi2 = 10.;
      double ndof = 10.;
      int Charge =1;
      v3.SetPtEtaPhi(pt, direction.eta() , direction.phi());
      const reco::TrackBase::Vector v(v3.x(), v3.y(), v3.z());
      const reco::TrackBase::Point vertex1(vertex.x(),vertex.y(), vertex.z());
      reco::Track psuedo_trk(chi2,ndof,vertex1,v , Charge,reco::Track::CovarianceMatrix());
      reco::IsoDeposit depTrk = muIsoExtractorTrack_->deposit(iEvent, iSetup, psuedo_trk );
      std::vector<reco::IsoDeposit> caloDeps = muIsoExtractorCalo_->deposits(iEvent, iSetup, psuedo_trk);
      //reco::IsoDeposit depJet = muIsoExtractorJet_->deposit(iEvent, iSetup, psuedo_trk);
      reco::IsoDeposit depEcal = caloDeps.at(0); 
      reco::IsoDeposit depHcal = caloDeps.at(1);
      reco::IsoDeposit depHo   = caloDeps.at(2);
      // std::cout <<"ecal  "<<depEcal.depositWithin(0.3)<<std::endl;
      
      ran_ecalIso03_mu  ->push_back(depEcal.depositWithin(0.3)  );
      ran_hcalIso03_mu  ->push_back(depHcal.depositWithin(0.3)  );
      ran_hoIso03_mu    ->push_back(depHo.depositWithin(0.3)    );
      ran_trkIso03_mu   ->push_back(depTrk.depositWithin(0.3)   );
      
      LorentzVector trks_p4    (directionAtVertex.x(),directionAtVertex.y(),directionAtVertex.z(), 0   );
      LorentzVector trksecal_p4(superCluster.x()     ,superCluster.y(),     superCluster.z(),      0   );
      ran_trksp4            ->push_back(trks_p4);
      ran_trksecalp4        ->push_back(trksecal_p4);
    }
  
  
  
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(ran_ecalIso03_egamma             ,branchprefix+"ecalIso03egamma"             );
  iEvent.put(ran_hcalIso03_egamma             ,branchprefix+"hcalIso03egamma"             );
  iEvent.put(ran_hcalD1Iso03_egamma           ,branchprefix+"hcalD1Iso03egamma"           );
  iEvent.put(ran_hcalD2Iso03_egamma           ,branchprefix+"hcalD2Iso03egamma"           );
  iEvent.put(ran_trkIso03_egamma              ,branchprefix+"trkIso03egamma"              );

  iEvent.put(ran_ecalIso03_mu                 ,branchprefix+"ecalIso03mu"                 );
  iEvent.put(ran_hcalIso03_mu                 ,branchprefix+"hcalIso03mu"                 );
  iEvent.put(ran_hoIso03_mu                   ,branchprefix+"hoIso03mu"                   );
  iEvent.put(ran_trkIso03_mu                  ,branchprefix+"trkIso03mu"                  );

  iEvent.put(ran_srFlag                       ,branchprefix+"srflag"                      );
  iEvent.put(ran_towerEmEt                    ,branchprefix+"toweremet"                   );
  iEvent.put(ran_towerHadEt                   ,branchprefix+"towerhadet"                  );
  iEvent.put(ran_dRClosestTower               ,branchprefix+"drclosesttower"              );
  iEvent.put(ran_dRClosestTowerEmEt           ,branchprefix+"drclosesttoweremet"          );

  iEvent.put(ran_trksp4                       ,branchprefix+"trksp4"                      );
  iEvent.put(ran_trksecalp4                   ,branchprefix+"trksecalp4"                 );
 

}
double RandomConeIsoMaker::signedRnd(double number, const double &min, const double &max)
{

	if (number > 0.5) {
		number = (number - 0.5) * 2.0;
		return (number * (max - min)) + min;
	}

	else {
		number = number * -2.0;
		return (number * (max - min)) - min;
	}

}

std::pair<int, double> RandomConeIsoMaker::getTrackIso(const reco::TrackBase::Point &beamspot,
		const reco::TrackCollection *trackCollection, 
		const math::XYZVector &tmpElectronMomentumAtVtx) const  
{
	int counter  = 0 ;
	double ptSum = 0.;

	for ( reco::TrackCollection::const_iterator itrTr  = (*trackCollection).begin() ; 
			itrTr != (*trackCollection).end()   ; 
			++itrTr ) {
		math::XYZVector tmpTrackMomentumAtVtx = (*itrTr).momentum () ; 

		double this_pt  = (*itrTr).pt();
		if ( this_pt < ptMin_ ) 
			continue;

		// note that pseudo directions are with respect to a 0, 0, 0 vertex
		// if (fabs( (*itrTr).dz() - (*tmpTrack).dz() ) > maxVtxDist_ )
		if (fabs( (*itrTr).dz() - 0.0 ) > maxVtxDist_ )
			continue;
		if (fabs( (*itrTr).dxy(beamspot) ) > drb_   )
			continue;
		double dr = ROOT::Math::VectorUtil::DeltaR(tmpTrackMomentumAtVtx,tmpElectronMomentumAtVtx) ;
		if ( fabs(dr) < extRadius_ && 
				fabs(dr) >= intRadius_ )
		{
			++counter ;
			ptSum += this_pt;
		}
	}//end loop over tracks                 

	std::pair<int,double> retval;
	retval.first  = counter;
	retval.second = ptSum;

	return retval;
}
//define this as a plug-in
DEFINE_FWK_MODULE(RandomConeIsoMaker);
