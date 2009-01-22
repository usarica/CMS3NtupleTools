// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      SCMaker
// 
/**\class SCMaker SCMaker.cc CMS2/NtupleMaker/src/SCMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
#ifndef CMS2_SCMAKER_H
#define CMS2_SCMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

//
// class declaration
//

class SCMaker : public edm::EDProducer {
public:
     explicit SCMaker (const edm::ParameterSet&);
  
private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

     math::XYZTLorentzVector initP4(const math::XYZPoint &pvPos,                 
                                    const reco::SuperCluster &sc);

     // ----------member data ---------------------------

     // preselection cuts
     double scEtMin_;

     // supercluster input collections
     edm::InputTag scInputTag_EE_;
     edm::InputTag scInputTag_EB_;
     std::vector<edm::InputTag> scInputTags_;

     // rechit input collections
     edm::InputTag hcalRecHitsInputTag_HBHE_;
     edm::InputTag ecalRecHitsInputTag_EE_;
     edm::InputTag ecalRecHitsInputTag_EB_;

     // primary vertex collection
     edm::InputTag primaryVertexInputTag_;

     // access to geometry
     unsigned long long cachedCaloGeometryID_;
     edm::ESHandle<CaloGeometry> caloGeometry_;

};

#endif

