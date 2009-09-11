// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      CaloTowerMaker
// 
/**\class CaloTowerMaker CaloTowerMaker.cc CMS2/NtupleMaker/src/CaloTowerMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
#ifndef CMS2_CALOTOWERMAKER_H
#define CMS2_CALOTOWERMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
// #include "RecoEgamma/EgammaTools/interface/HoECalculator.h"
// #include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
// #include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

//
// class declaration
//

class CaloTowerMaker : public edm::EDProducer {
public:
     explicit CaloTowerMaker (const edm::ParameterSet&);
  
private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

//      math::XYZTLorentzVector initP4(const math::XYZPoint &pvPos,                 
//                                     const reco::SuperCluster &sc);

     // ----------member data ---------------------------

//      // preselection cuts
//      double scEtMin_;

     // primary vertex collection
     edm::InputTag primaryVertexInputTag_;

     // Calo Tower collection
     edm::InputTag caloTowersInputTag_;

//      // access to geometry
//      unsigned long long cachedCaloGeometryID_;
//      edm::ESHandle<CaloGeometry> caloGeometry_;

};

#endif

