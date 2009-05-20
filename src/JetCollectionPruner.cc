// -*- C++ -*-
//
// Package:    JetCollectionPruner
// Class:      JetCollectionPruner
// 
/**\class JetCollectionPruner JetCollectionPruner.cc CMS2/JetCollectionPruner/src/JetCollectionPruner.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "CMS2/NtupleMaker/interface/JetCollectionPruner.h"


bool sortJetsByPt(reco::CaloJet jet1,
		  reco::CaloJet jet2) {
  return jet1.pt() > jet2.pt();
}


//
// constructors and destructor
//
JetCollectionPruner::JetCollectionPruner(const edm::ParameterSet& iConfig)
{

     // define what to produce
  produces<reco::CaloJetCollection>();

  // get the input collection of electrons
  inputUncorrectedJetCollection_ = iConfig.getParameter<edm::InputTag>("inputUncorrectedJetCollection");
  uncorrectedJetPtCut_ = iConfig.getParameter<double>("uncorrectedJetPtCut");

}

void JetCollectionPruner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

     // get the primary vertices
     edm::Handle<reco::CaloJetCollection> jetsHandle;
     try {
        iEvent.getByLabel(inputUncorrectedJetCollection_, jetsHandle);
     }
     catch ( cms::Exception& ex ) {
        edm::LogError("JetCollectionPrunerError") << "Error! can't get the input recoJets";
     }
     // new collection to put in the event after cleaning
     std::auto_ptr<reco::CaloJetCollection> outCollection (new reco::CaloJetCollection);

     for(reco::CaloJetCollection::const_iterator jet_it = jetsHandle->begin();
	 jet_it != jetsHandle->end(); jet_it++) {

       if(jet_it->pt() > uncorrectedJetPtCut_)
	 outCollection->push_back(*jet_it);

     }

     //sort the jets using the above predicate.
     //they should be sorted, but being paranoid
     std::sort(outCollection->begin(), outCollection->end(), sortJetsByPt);

     // put the cleaned collection in the event
     iEvent.put(outCollection);

}

// ------------ method called once each job just before starting event loop  ------------
void 
JetCollectionPruner::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetCollectionPruner::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetCollectionPruner);

