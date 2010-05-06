// -*- C++ -*-
//
// Package:    TrackJetCollectionPruner
// Class:      TrackJetCollectionPruner
// 
/**\class TrackJetCollectionPruner TrackJetCollectionPruner.cc CMS2/TrackJetCollectionPruner/src/TrackJetCollectionPruner.cc

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
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "CMS2/NtupleMaker/interface/TrackJetCollectionPruner.h"

bool sortJetsByPt(reco::TrackJet jet1, reco::TrackJet jet2) {
  return jet1.pt() > jet2.pt();
}

//
// constructors and destructor
//
TrackJetCollectionPruner::TrackJetCollectionPruner(const edm::ParameterSet& iConfig)
{
  // define what to produce
  produces<reco::TrackJetCollection>();

  // get the input collection of electrons
  inputUncorrectedJetCollection_ = iConfig.getParameter<edm::InputTag>("inputUncorrectedJetCollection");
  uncorrectedJetPtCut_           = iConfig.getParameter<double>       ("uncorrectedJetPtCut"          );

}

void TrackJetCollectionPruner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get the primary vertices
  edm::Handle<reco::TrackJetCollection> jetsHandle;
  try {
    iEvent.getByLabel(inputUncorrectedJetCollection_, jetsHandle);
  }
  catch ( cms::Exception& ex ) {
    edm::LogError("TrackJetCollectionPrunerError") << "Error! can't get the input recoJets";
  }

  // new collection to put in the event after cleaning
  std::auto_ptr<reco::TrackJetCollection> outCollection (new reco::TrackJetCollection);

  for(reco::TrackJetCollection::const_iterator jet_it = jetsHandle->begin(); jet_it != jetsHandle->end(); jet_it++) {

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
void TrackJetCollectionPruner::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void TrackJetCollectionPruner::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackJetCollectionPruner);

