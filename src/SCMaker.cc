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

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;

//
// class decleration
//

//
// constructors and destructor
//
SCMaker::SCMaker(const edm::ParameterSet& iConfig)
{
     produces<std::vector<float> >("scsenergy").setBranchAlias("scs_energy"); // sc energy

     scInputTag_EE_ = iConfig.getParameter<edm::InputTag>("scInputTag_EE");
     scInputTag_EB_ = iConfig.getParameter<edm::InputTag>("scInputTag_EB");

}

void SCMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;

     // make vectors to hold the information
     std::auto_ptr<std::vector<float> > vector_scs_energy (new std::vector<float>);      

     // get superclusters

     Handle<reco::SuperClusterCollection> scHandle_EE;
     try { 
        iEvent.getByLabel(scInputTag_EE_, scHandle_EE);
     } catch ( cms::Exception& ex ) {
        edm::LogError("SCMakerError") 
           << "Error! can't get the SuperClusters " 
           << scInputTag_EE_.label() ;
     }    

     const reco::SuperClusterCollection *scCollection_EE = scHandle_EE.product();

     for (reco::SuperClusterCollection::const_iterator sc = scCollection_EE->begin(); 
          				sc != scCollection_EE->end(); ++sc) {
	  // fill vectorsw
	  vector_scs_energy->push_back( sc->energy() );
     }

     // store vectors
     iEvent.put(vector_scs_energy, "scsenergy");

}

// ------------ method called once each job just before starting event loop  ------------
void 
SCMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(SCMaker);

