// -*- C++ -*-
//
// Package:    PFClusterMaker
// Class:      PFClusterMaker
// 
/**\class PFClusterMaker PFClusterMaker.cc CMS2/NtupleMaker/src/PFClusterMaker.cc

   Description: fill collection of PFClusters

   Implementation:
   - extract and fill variables
*/
//
// Original Ben Hooberman
// Created:  Wed Mar  24 12:23:38 CDT 2010
// 
//
//


// system include files
#include <memory>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "CMS2/NtupleMaker/interface/PFClusterMaker.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "TMath.h"

typedef ROOT::Math::PositionVector3D<ROOT::Math::CylindricalEta3D<Double32_t> > REPPoint;
//
// class declaration
//

//
// constructors and destructor
//
using namespace edm;
using namespace std;

PFClusterMaker::PFClusterMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // product of this EDProducer
  produces<float>	       (branchprefix+"met").     setBranchAlias(aliasprefix_+"_met");
  produces<float>	       (branchprefix+"sumet").   setBranchAlias(aliasprefix_+"_sumet");
  produces<float>	       (branchprefix+"metphi").  setBranchAlias(aliasprefix_+"_metphi");
  produces<std::vector<float> >(branchprefix+"energy").  setBranchAlias(aliasprefix_+"_energy");
  produces<std::vector<float> >(branchprefix+"et").      setBranchAlias(aliasprefix_+"_et");
  produces<std::vector<float> >(branchprefix+"eta").     setBranchAlias(aliasprefix_+"_eta");
  produces<std::vector<float> >(branchprefix+"phi").     setBranchAlias(aliasprefix_+"_phi");
  produces<std::vector<int> >  (branchprefix+"layer").   setBranchAlias(aliasprefix_+"_layer");
  
  // parameters from configuration
  inputTagPFClustersECAL_  = iConfig.getParameter<edm::InputTag>("PFClustersECAL");
  inputTagPFClustersHCAL_  = iConfig.getParameter<edm::InputTag>("PFClustersHCAL");
  inputTagPFClustersHFEM_  = iConfig.getParameter<edm::InputTag>("PFClustersHFEM");
  inputTagPFClustersHFHAD_ = iConfig.getParameter<edm::InputTag>("PFClustersHFHAD");


}

PFClusterMaker::~PFClusterMaker() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void PFClusterMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {


  std::auto_ptr<float> pfc_met		        (new float);
  std::auto_ptr<float> pfc_sumet		(new float);
  std::auto_ptr<float> pfc_metphi		(new float);
  std::auto_ptr<std::vector<float> > pfc_energy (new std::vector<float>);
  std::auto_ptr<std::vector<float> > pfc_et     (new std::vector<float>);
  std::auto_ptr<std::vector<float> > pfc_eta    (new std::vector<float>);
  std::auto_ptr<std::vector<float> > pfc_phi    (new std::vector<float>);
  std::auto_ptr<std::vector<int> >   pfc_layer  (new std::vector<int>);

  bool found;

  //ECAL PFClusters
  edm::Handle< reco::PFClusterCollection > clustersECAL;
  found = iEvent.getByLabel(inputTagPFClustersECAL_, clustersECAL);      

  if(!found )
    LogError("PFBlockProducer")<<" cannot get ECAL clusters: "
                               <<inputTagPFClustersECAL_<<endl;

  //HCAL PFClusters
  edm::Handle< reco::PFClusterCollection > clustersHCAL;
  found = iEvent.getByLabel(inputTagPFClustersHCAL_, clustersHCAL);
  
  if(!found )
    LogError("PFBlockProducer")<<" cannot get HCAL clusters: "
                               <<inputTagPFClustersHCAL_<<endl;

  //HFEM PFClusters
  edm::Handle< reco::PFClusterCollection > clustersHFEM;
  found = iEvent.getByLabel(inputTagPFClustersHFEM_, clustersHFEM);      
  
  if(!found )
    LogError("PFBlockProducer")<<" cannot get HFEM clusters: "
                               <<inputTagPFClustersHFEM_<<endl;
  
  //HFHAD PFClusters
  edm::Handle< reco::PFClusterCollection > clustersHFHAD;
  found = iEvent.getByLabel(inputTagPFClustersHFHAD_, clustersHFHAD);      
  
  if(!found )
    LogError("PFBlockProducer")<<" cannot get HFHAD clusters: "
                               <<inputTagPFClustersHFHAD_<<endl;

  float met_x = 0.;
  float met_y = 0.;
  float sumet = 0.;

  //---------------ECAL PFClusters---------------
  
  for (reco::PFClusterCollection::const_iterator it = clustersECAL->begin(); it != clustersECAL->end(); it++){

    const math::XYZPoint&  cluster_pos = it->position();
    double et = it->energy() / cosh( cluster_pos.eta() ); 

    pfc_energy	->push_back( it->energy()	);
    pfc_et	->push_back( et			);
    pfc_eta	->push_back( cluster_pos.eta()	);
    pfc_phi	->push_back( cluster_pos.phi()	);
    pfc_layer	->push_back( it->layer()        );

    met_x    -= et * cos(cluster_pos.phi());
    met_y    -= et * sin(cluster_pos.phi());
    sumet    += et;      
  }

  //---------------HCAL PFClusters---------------
  
  for (reco::PFClusterCollection::const_iterator it = clustersHCAL->begin(); it != clustersHCAL->end(); it++){

    if ( it->layer() == PFLayer::HCAL_BARREL2) continue; //skip HO

    const math::XYZPoint&  cluster_pos = it->position();
    double et = it->energy() / cosh( cluster_pos.eta() ); 

    pfc_energy	->push_back( it->energy()	);
    pfc_et	->push_back( et			);
    pfc_eta	->push_back( cluster_pos.eta()	);
    pfc_phi	->push_back( cluster_pos.phi()	);
    pfc_layer	->push_back( it->layer()        );
 
    met_x    -= et * cos(cluster_pos.phi());
    met_y    -= et * sin(cluster_pos.phi());
    sumet    += et;   
    
  }

  //---------------HFEM PFClusters---------------
  
  for (reco::PFClusterCollection::const_iterator it = clustersHFEM->begin(); it != clustersHFEM->end(); it++){

    const math::XYZPoint&  cluster_pos = it->position();
    double et = it->energy() / cosh( cluster_pos.eta() ); 

    pfc_energy	->push_back( it->energy()	);
    pfc_et	->push_back( et			);
    pfc_eta	->push_back( cluster_pos.eta()	);
    pfc_phi	->push_back( cluster_pos.phi()	);
    pfc_layer	->push_back( it->layer()        );

    met_x    -= et * cos(cluster_pos.phi());
    met_y    -= et * sin(cluster_pos.phi());
    sumet    += et;   
  }

  //---------------HFHAD PFClusters---------------
  
  for (reco::PFClusterCollection::const_iterator it = clustersHFHAD->begin(); it != clustersHFHAD->end(); it++){

    const math::XYZPoint&  cluster_pos = it->position();
    double et = it->energy() / cosh( cluster_pos.eta() ); 

    pfc_energy	->push_back( it->energy()	);
    pfc_et	->push_back( et			);
    pfc_eta	->push_back( cluster_pos.eta()	);
    pfc_phi	->push_back( cluster_pos.phi()	);
    pfc_layer	->push_back( it->layer()        );

    met_x    -= et * cos(cluster_pos.phi());
    met_y    -= et * sin(cluster_pos.phi());
    sumet    += et;   
  }

  *pfc_met    = sqrt( met_x * met_x + met_y * met_y );
  *pfc_sumet  = sumet;
  *pfc_metphi = atan2( met_y , met_x );

  //---------------put containers into event---------------

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(pfc_met,      branchprefix+"met");
  iEvent.put(pfc_sumet,    branchprefix+"sumet");
  iEvent.put(pfc_metphi,   branchprefix+"metphi");
  iEvent.put(pfc_energy,   branchprefix+"energy");
  iEvent.put(pfc_et,       branchprefix+"et");
  iEvent.put(pfc_eta,      branchprefix+"eta");
  iEvent.put(pfc_phi,      branchprefix+"phi");
  iEvent.put(pfc_layer,    branchprefix+"layer");
  
}

// ------------ method called once each job just before starting event loop  ------------
void PFClusterMaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void PFClusterMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFClusterMaker);
