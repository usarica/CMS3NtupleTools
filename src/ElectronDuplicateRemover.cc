// -*- C++ -*-
//
// Package:    ElectronDuplicateRemover
// Class:      ElectronDuplicateRemover
// 
/**\class ElectronDuplicateRemover ElectronDuplicateRemover.cc CMS2/ElectronDuplicateRemover/src/ElectronDuplicateRemover.cc

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

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "CMS2/NtupleMaker/interface/ElectronDuplicateRemover.h"

//
// constructors and destructor
//
ElectronDuplicateRemover::ElectronDuplicateRemover(const edm::ParameterSet& iConfig)
{

     // define what to produce
     produces<reco::GsfElectronCollection>();

     // get the input collection of electrons
     electronsInputTag_ = iConfig.getParameter<edm::InputTag>("electronsInputTag");

}

void ElectronDuplicateRemover::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

     // get the primary vertices
     edm::Handle<reco::GsfElectronCollection> electronHandle;
     try {
        iEvent.getByLabel(electronsInputTag_, electronHandle);
     }
     catch ( cms::Exception& ex ) {
        edm::LogError("ElectronDuplicateRemoverError") << "Error! can't get the input electrons";
     }
     const reco::GsfElectronCollection *elCollection = electronHandle.product();
     reco::GsfElectronCollection sortedCollection(elCollection->size());
     std::copy(elCollection->begin(), elCollection->end(), sortedCollection.begin());

     // new collection to put in the event after cleaning
     std::auto_ptr<reco::GsfElectronCollection> cleanedCollection (new reco::GsfElectronCollection);

     reco::GsfElectronCollection::iterator e1, e2;
     std::sort(sortedCollection.begin(), sortedCollection.end(), BetterElectron());

     // resolve when e/g SC is found
     for(e1 = sortedCollection.begin(); e1 != sortedCollection.end(); ++e1)
     {
          reco::SuperClusterRef scRef1 = e1->superCluster();
          for (e2 = e1, ++e2; e2 != sortedCollection.end();)
          {
               reco::SuperClusterRef scRef2 = e2->superCluster();
               if (scRef1 == scRef2 || e1->gsfTrack() == e2->gsfTrack())
               {
                    e2 = sortedCollection.erase(e2);
               }
               else
               { ++e2; }
          }
          cleanedCollection->push_back(*e1);
     }

     // put the cleaned collection in the event
     iEvent.put(cleanedCollection);

}

/*
bool ElectronDuplicateRemover::betterElectron(const reco::GsfElectron &e1, const reco::GsfElectron &e2)
{
     return (fabs(e1.eSuperClusterOverP() - 1) < fabs(e2.eSuperClusterOverP() - 1));
}
*/

// ------------ method called once each job just before starting event loop  ------------
void 
ElectronDuplicateRemover::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronDuplicateRemover::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronDuplicateRemover);

