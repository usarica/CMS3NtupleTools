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

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "CMS2/NtupleMaker/interface/CaloTowerMaker.h"

//#include "DataFormats/CaloRecHit/interface/CaloID.h"
//#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;

//
// class decleration
//

//
// constructors and destructor
//
CaloTowerMaker::CaloTowerMaker(const edm::ParameterSet& iConfig)
{
  //-----------------------------------------------------------------------------------------
  //--- which of these will we need? -- uncommment what will be dumped ----------------------
  //-----------------------------------------------------------------------------------------
//OUT//   // getters
//OUT//   CaloTowerDetId id() const { return id_; }
//OUT//   const std::vector<DetId>& constituents() const { return constituents_; }
//OUT//   size_t constituentsSize() const { return constituents_.size(); }
//OUT//   DetId constituent( size_t i ) const { return constituents_[ i ]; }
//OUT
//OUT//   // energy contributions from different detectors
//OUT//   // energy in HO ("outerEnergy")is not included in "hadEnergy"
//OUT   double emEnergy() const { return emE_ ; }
//OUT   double hadEnergy() const { return hadE_ ; }
//OUT   double outerEnergy() const { return (id_.ietaAbs()<16)? outerE_ : 0.0; }
//OUT
//OUT//   // transverse energies wrt to vtx (0,0,0)
//OUT   double emEt() const { return emE_ * sin( theta() ); }
//OUT   double hadEt() const { return hadE_ * sin( theta() ); }
//OUT   double outerEt() const { return (id_.ietaAbs()<16)? outerE_ * sin( theta() ) : 0.0; }
//OUT
//OUT
//OUT//   // preserve the inherited default accessors where applicable
//OUT//   // (user gets default p4 wrt to vtx (0,0,0) using p4(), etc.
//OUT//   using LeafCandidate::p4;
//OUT//   using LeafCandidate::p;
//OUT//   using LeafCandidate::et; 
//OUT
//OUT//   // recalculated wrt user provided vertex Z position;
//OUT//    math::PtEtaPhiMLorentzVector p4(double vtxZ) const;
//OUT//    double p (double vtxZ) const { return p4(vtxZ).P(); }
//OUT//    double et(double vtxZ) const { return p4(vtxZ).Et(); }
//OUT
//OUT//    double emEt(double vtxZ)  const { return  emE_ * sin(p4(vtxZ).theta()); }
//OUT//    double hadEt(double vtxZ) const { return  hadE_ * sin(p4(vtxZ).theta()); }
//OUT//    double outerEt(double vtxZ) const { return (id_.ietaAbs()<16)? outerE_ * sin(p4(vtxZ).theta()) : 0.0; }
//OUT
//OUT//   // recalculated wrt vertex provided as 3D point
//OUT   math::PtEtaPhiMLorentzVector p4(Point v) const;
//OUT   double p (Point v) const { return p4(v).P(); }
//OUT   double et(Point v) const { return p4(v).Et(); }
//OUT   double emEt(Point v)  const { return  emE_ * sin(p4(v).theta()); }
//OUT   double hadEt(Point v) const { return  hadE_ * sin(p4(v).theta()); }
//OUT   double outerEt(Point v) const { return (id_.ietaAbs()<16)? outerE_ * sin(p4(v).theta()) : 0.0; }
//OUT
//OUT//   // Access to p4 comming from HO alone: requested by JetMET to add/subtract HO contributions
//OUT//   // to the tower for cases when the tower collection was created without/with HO   
//OUT//   math::PtEtaPhiMLorentzVector p4_HO() const;  
//OUT//   math::PtEtaPhiMLorentzVector p4_HO(double vtxZ) const;
//OUT//   math::PtEtaPhiMLorentzVector p4_HO(Point v) const;
//OUT
//OUT//   // the reference poins in ECAL and HCAL for direction determination
//OUT//   // algorithm and parameters for selecting these points are set in the CaloTowersCreator
//OUT   const GlobalPoint& emPosition()  const { return emPosition_ ; }
//OUT   const GlobalPoint& hadPosition() const { return hadPosition_ ; }
//OUT
//OUT//   int emLvl1() const { return emLvl1_; }
//OUT//   int hadLv11() const { return hadLvl1_; }
//OUT
//OUT//   // energy contained in depths>1 in the HE for 18<|iEta|<29
//OUT//   double hadEnergyHeOuterLayer() const { return (id_.ietaAbs()<18 || id_.ietaAbs()>29)? 0 : outerE_; }
//OUT//   double hadEnergyHeInnerLayer() const { return (id_.ietaAbs()<18 || id_.ietaAbs()>29)? 0 : hadE_ - outerE_; }
//OUT
//OUT//   // time (ns) in ECAL/HCAL components of the tower based on weigted sum of the times in the contributing RecHits
//OUT   float ecalTime() const { return float(ecalTime_) * 0.01; }
//OUT   float hcalTime() const { return float(hcalTime_) * 0.01; }
//OUT
//OUT//   // position information on the tower
//OUT//   int ieta() const { return id_.ieta(); }
//OUT//   int ietaAbs() const { return id_.ietaAbs(); }
//OUT//   int iphi() const { return id_.iphi(); }
//OUT//   int zside() const { return id_.zside(); }
//OUT
//OUT   int numCrystals() const; 
//OUT
//OUT//   // methods to retrieve status information from the CaloTower:
//OUT//   // number of bad/recovered/problematic cells in the tower
//OUT//   // separately for ECAL and HCAL
//OUT
//OUT   uint numBadEcalCells() const { return (twrStatusWord_ & 0x1F); }
//OUT   uint numRecoveredEcalCells() const { return ((twrStatusWord_ >> 5) & 0x1F); }
//OUT   uint numProblematicEcalCells() const { return ((twrStatusWord_ >> 10) & 0x1F); }
//OUT
//OUT   uint numBadHcalCells() const { return ( (twrStatusWord_ >> 15)& 0x7); }
//OUT   uint numRecoveredHcalCells() const { return ((twrStatusWord_ >> 18) & 0x7); }
//OUT   uint numProblematicHcalCells() const { return ((twrStatusWord_ >> 21) & 0x7); }
//OUT
//OUT//   // the status word itself
//OUT//   uint32_t towerStatusWord() const { return twrStatusWord_; }
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  // number of towers in the event
  produces<unsigned int>("evtntwrs").setBranchAlias("evt_ntwrs");

  produces<std::vector<float> >("twrseta").setBranchAlias("twrs_eta");
  produces<std::vector<float> >("twrsphi").setBranchAlias("twrs_phi");

  // energy contributions from different detectors
  // energy in HO ("outerEnergy")is not included in "hadEnergy"
  //   double emEnergy() const { return emE_ ; }
  produces<std::vector<float> >("twrsemEnergy").setBranchAlias("twrs_emEnergy");
  //   double hadEnergy() const { return hadE_ ; }
  produces<std::vector<float> >("twrshadEnergy").setBranchAlias("twrs_hadEnergy");
  //   double outerEnergy() const { return (id_.ietaAbs()<16)? outerE_ : 0.0; }
  produces<std::vector<float> >("twrsouterEnergy").setBranchAlias("twrs_outerEnergy");

   // transverse energies wrt to vtx (0,0,0)
  //   double emEt() const { return emE_ * sin( theta() ); }
  produces<std::vector<float> >("twrsemEt").setBranchAlias("twrs_emEt");
  //   double hadEt() const { return hadE_ * sin( theta() ); }
  produces<std::vector<float> >("twrshadEt").setBranchAlias("twrs_hadEt");
  //   double outerEt() const { return (id_.ietaAbs()<16)? outerE_ * sin( theta() ) : 0.0; }
  produces<std::vector<float> >("twrsouterEt").setBranchAlias("twrs_outerEt");

   // recalculated wrt vertex provided as 3D point
  //   math::PtEtaPhiMLorentzVector p4(Point v) const;
  //   double p (Point v) const { return p4(v).P(); }
  produces<std::vector<float> >("twrspcorr").setBranchAlias("twrs_pcorr");
  //   double et(Point v) const { return p4(v).Et(); }
  produces<std::vector<float> >("twrsetcorr").setBranchAlias("twrs_etcorr");
  //   double emEt(Point v)  const { return  emE_ * sin(p4(v).theta()); }
  produces<std::vector<float> >("twrsemEtcorr").setBranchAlias("twrs_emEtcorr");
  //   double hadEt(Point v) const { return  hadE_ * sin(p4(v).theta()); }
  produces<std::vector<float> >("twrshadEtcorr").setBranchAlias("twrs_hadEtcorr");
  //   double outerEt(Point v) const { return (id_.ietaAbs()<16)? outerE_ * sin(p4(v).theta()) : 0.0; }
  produces<std::vector<float> >("twrsouterEtcorr").setBranchAlias("twrs_outerEtcorr");

//    // the reference poins in ECAL and HCAL for direction determination
//    // algorithm and parameters for selecting these points are set in the CaloTowersCreator
//    const GlobalPoint& emPosition()  const { return emPosition_ ; }
//    const GlobalPoint& hadPosition() const { return hadPosition_ ; }

   // time (ns) in ECAL/HCAL components of the tower based on weigted sum of the times in the contributing RecHits
  //   float ecalTime() const { return float(ecalTime_) * 0.01; }
  produces<std::vector<float> >("twrsecalTime").setBranchAlias("twrs_ecalTime");
  //   float hcalTime() const { return float(hcalTime_) * 0.01; }
  produces<std::vector<float> >("twrshcalTime").setBranchAlias("twrs_hcalTime");

  // methods to retrieve status information from the CaloTower:
  // number of bad/recovered/problematic cells in the tower
  // separately for ECAL and HCAL
  //  uint numBadEcalCells() const { return (twrStatusWord_ & 0x1F); }
  produces<std::vector<unsigned int> >("twrsnumBadEcalCells").setBranchAlias("twrs_numBadEcalCells");
  //  uint numRecoveredEcalCells() const { return ((twrStatusWord_ >> 5) & 0x1F); }
  produces<std::vector<unsigned int> >("twrsnumRecoveredEcalCells").setBranchAlias("twrs_numRecoveredEcalCells");
  //  uint numProblematicEcalCells() const { return ((twrStatusWord_ >> 10) & 0x1F); }
  produces<std::vector<unsigned int> >("twrsnumProblematicEcalCells").setBranchAlias("twrs_numProblematicEcalCells");

  //  uint numBadHcalCells() const { return ( (twrStatusWord_ >> 15)& 0x7); }
  produces<std::vector<unsigned int> >("twrsnumBadHcalCells").setBranchAlias("twrs_numBadHcalCells");
  //  uint numRecoveredHcalCells() const { return ((twrStatusWord_ >> 18) & 0x7); }
  produces<std::vector<unsigned int> >("twrsnumRecoveredHcalCells").setBranchAlias("twrs_numRecoveredHcalCells");
  //  uint numProblematicHcalCells() const { return ((twrStatusWord_ >> 21) & 0x7); }
  produces<std::vector<unsigned int> >("twrsnumProblematicHcalCells").setBranchAlias("twrs_numProblematicHcalCells");

  // the energy of the max energy crystal in the tower
  produces<std::vector<float> >("twrsemMax").setBranchAlias("twrs_emMax");
  // the energy in 3x3 crystals centred on the max energy crystal
  produces<std::vector<float> >("twrsem3x3").setBranchAlias("twrs_em3x3");
  // as above for 5x5 crystals
  produces<std::vector<float> >("twrsem5x5").setBranchAlias("twrs_em5x5");
  
  // add superclusters to the ntuple if they have ET > scEtMin_
  //   scEtMin_ = iConfig.getParameter<double>("scEtMin");

   // input Tags
   primaryVertexInputTag_ = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
   caloTowersInputTag_ = iConfig.getParameter<edm::InputTag>("caloTowersInputTag");
   ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
   ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");

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

  // ecal cluster shape variables
  // do not use the lazy tools because need to get the hits anyway
  EcalClusterTools clusterTools;

//   // get hoe variable
//   HoECalculator hoeCalc(caloGeometry_);

  std::auto_ptr<unsigned int> evt_ntwrs (new unsigned int);

  std::auto_ptr<std::vector<float> > vector_twrs_eta (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_twrs_phi (new std::vector<float>);

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

  std::auto_ptr<std::vector<float> > vector_twrs_ecalTime (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_twrs_hcalTime (new std::vector<float>);

  std::auto_ptr<std::vector<unsigned int> > vector_twrs_numBadEcalCells (new std::vector<unsigned int>);
  std::auto_ptr<std::vector<unsigned int> > vector_twrs_numRecoveredEcalCells (new std::vector<unsigned int>);
  std::auto_ptr<std::vector<unsigned int> > vector_twrs_numProblematicEcalCells (new std::vector<unsigned int>);

  std::auto_ptr<std::vector<unsigned int> > vector_twrs_numBadHcalCells (new std::vector<unsigned int>);
  std::auto_ptr<std::vector<unsigned int> > vector_twrs_numRecoveredHcalCells (new std::vector<unsigned int>);
  std::auto_ptr<std::vector<unsigned int> > vector_twrs_numProblematicHcalCells (new std::vector<unsigned int>);

  std::auto_ptr<std::vector<float> > vector_twrs_emMax (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_twrs_em3x3 (new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_twrs_em5x5 (new std::vector<float>);

  *evt_ntwrs = 0;

  for (CaloTowerCollection::const_iterator j = calotower->begin();j !=calotower->end(); j++) {

    *evt_ntwrs +=1;
    
//     std::cout << *j << std::endl;
//     std::cout << "ENERGY HAD " << j->hadEnergy()<< " ENERGY EM " <<j->emEnergy() 
//          << " ETA " <<j->eta() << " PHI " <<j->phi() << std::endl;

    vector_twrs_eta->push_back(j->eta());
    vector_twrs_phi->push_back(j->phi());

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
    float emMax = 0.0;
    float emE = 0.0;
    DetId emMaxId(0);

    // loop on detids in the tower
    for (size_t i = 0; i < towerDetIds.size(); ++i) {

      // find the energy of this detId if it is in the ecal
      if (towerDetIds[i].subdetId() == EcalEndcap)
         emE = clusterTools.recHitEnergy(towerDetIds[i], recHitsEE);
      if (towerDetIds[i].subdetId() == EcalBarrel) 
         emE = clusterTools.recHitEnergy(towerDetIds[i], recHitsEB);

      // compare with previous highest energy
      if (emE > emMax) {
           emMax = emE; 
	   emMaxId = towerDetIds[i];
        }
    }

    // get the e3x3 and e5x5 centred on this detid
    // - this will include crystals not in the tower if necessary
    float em3x3;
    float em5x5;
    // the function to compute the matrix energy takes a basic cluster as an argument
    // but uses nothing from it!?!
    reco::BasicCluster dummyCluster;
    if (emMaxId.subdetId() == EcalEndcap) {
        em3x3 = clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, emMaxId, -1, 1, -1, 1);
        em5x5 = clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, emMaxId, -2, 2, -2, 2);
    }
    if (emMaxId.subdetId() == EcalBarrel) {
        em3x3 = clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, emMaxId, -1, 1, -1, 1);
        em5x5 = clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, emMaxId, -2, 2, -2, 2);
    }

    vector_twrs_emMax->push_back(emMax);
    vector_twrs_em3x3->push_back(em3x3);
    vector_twrs_em5x5->push_back(em5x5);

  }

// 	vector_scs_sigmaEtaEta->push_back( sqrt(covariances[0]) );

  // put results into the event
  iEvent.put(evt_ntwrs, "evtntwrs");
  iEvent.put(vector_twrs_eta, "twrseta");
  iEvent.put(vector_twrs_phi, "twrsphi");
  iEvent.put(vector_twrs_emEnergy, "twrsemEnergy");
  iEvent.put(vector_twrs_hadEnergy, "twrshadEnergy");
  iEvent.put(vector_twrs_outerEnergy, "twrsouterEnergy");

  iEvent.put(vector_twrs_emEt, "twrsemEt");
  iEvent.put(vector_twrs_hadEt, "twrshadEt");
  iEvent.put(vector_twrs_outerEt, "twrsouterEt");

  iEvent.put(vector_twrs_pcorr, "twrspcorr");
  iEvent.put(vector_twrs_etcorr, "twrsetcorr");
  iEvent.put(vector_twrs_emEtcorr, "twrsemEtcorr");
  iEvent.put(vector_twrs_hadEtcorr, "twrshadEtcorr");
  iEvent.put(vector_twrs_outerEtcorr, "twrsouterEtcorr");

  iEvent.put(vector_twrs_ecalTime, "twrsecalTime");
  iEvent.put(vector_twrs_hcalTime, "twrshcalTime");

  iEvent.put(vector_twrs_numBadEcalCells, "twrsnumBadEcalCells");
  iEvent.put(vector_twrs_numRecoveredEcalCells, "twrsnumRecoveredEcalCells");
  iEvent.put(vector_twrs_numProblematicEcalCells, "twrsnumProblematicEcalCells");

  iEvent.put(vector_twrs_numBadHcalCells, "twrsnumBadHcalCells");
  iEvent.put(vector_twrs_numRecoveredHcalCells, "twrsnumRecoveredHcalCells");
  iEvent.put(vector_twrs_numProblematicHcalCells, "twrsnumProblematicHcalCells");

  iEvent.put(vector_twrs_emMax, "twrsemMax");
  iEvent.put(vector_twrs_em3x3, "twrsem3x3");
  iEvent.put(vector_twrs_em5x5, "twrsem5x5");

}

// math::XYZTLorentzVector CaloTowerMaker::initP4(const math::XYZPoint &pvPos, 
//                                         const reco::SuperCluster &sc)
// {

//   math::XYZVector scPos(sc.x(), sc.y(), sc.z());
//   math::XYZVector pvPosVec(pvPos.x(), pvPos.y(), pvPos.z());
//   math::XYZVector objPosition = scPos - pvPosVec;
//   double scale = sc.energy() / objPosition.R();
//   return math::XYZTLorentzVector(objPosition.x() * scale, 
// 				 objPosition.y() * scale, 
// 				 objPosition.z() * scale, 
// 				 sc.energy());
// }



// ------------ method called once each job just before starting event loop  ------------
void 
CaloTowerMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CaloTowerMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloTowerMaker);

