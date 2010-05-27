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
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "CMS2/NtupleMaker/interface/JetCollectionPruner.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

bool sortJetsByPt(reco::JPTJet jet1, reco::JPTJet jet2) {
     return jet1.pt() > jet2.pt();
}

//
// constructors and destructor
//
JetCollectionPruner::JetCollectionPruner(const edm::ParameterSet& iConfig)
{
     // define what to produce
     produces<reco::JPTJetCollection>("jpt");
     produces<reco::CaloJetCollection>("calojet");

     // get the input collection of electrons
     inputUncorrectedJetCollection_ = iConfig.getParameter<edm::InputTag>("inputUncorrectedJetCollection");
     uncorrectedJetPtCut_           = iConfig.getParameter<double>       ("uncorrectedJetPtCut"          );
}

void JetCollectionPruner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     // get the primary vertices
     edm::Handle<reco::JPTJetCollection> jetsHandle;
     try {
	  iEvent.getByLabel(inputUncorrectedJetCollection_, jetsHandle);
     }
     catch ( cms::Exception& ex ) {
	  edm::LogError("JetCollectionPrunerError") << "Error! can't get the input recoJets";
     }

     // new collection to put in the event after cleaning
     std::auto_ptr<reco::JPTJetCollection> outJPTcollection (new reco::JPTJetCollection);
     std::auto_ptr<reco::CaloJetCollection> outJETcollection (new reco::CaloJetCollection);

     for (reco::JPTJetCollection::const_iterator jet_it = jetsHandle->begin(); jet_it != jetsHandle->end(); jet_it++)
     {
	  if(jet_it->pt() > uncorrectedJetPtCut_)
	       outJPTcollection->push_back(*jet_it);
     }

     //sort the jets using the above predicate.
     //they should be sorted, but being paranoid
     std::sort(outJPTcollection->begin(), outJPTcollection->end(), sortJetsByPt);

     for (reco::JPTJetCollection::const_iterator jet_it = outJPTcollection->begin(); jet_it != outJPTcollection->end(); jet_it++)
     {
	  const reco::CaloJet& tmp_calojet = dynamic_cast<const reco::CaloJet&>(*(jet_it->getCaloJetRef()));
	  outJETcollection->push_back(tmp_calojet);
     }

     // put the cleaned collection in the event
     iEvent.put(outJPTcollection, "jpt");
     iEvent.put(outJETcollection, "calojet");
}

// ------------ method called once each job just before starting event loop  ------------
void JetCollectionPruner::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void JetCollectionPruner::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetCollectionPruner);

