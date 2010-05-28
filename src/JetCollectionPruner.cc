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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "CMS2/NtupleMaker/interface/JetCollectionPruner.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

bool sortJPTJetsByPt(reco::JPTJet jet1, reco::JPTJet jet2) {
  return jet1.pt() > jet2.pt();
}

bool sortPFJetsByPt(reco::PFJet jet1, reco::PFJet jet2) {
  return jet1.pt() > jet2.pt();
}

bool sortTrkJetsByPt(reco::TrackJet jet1, reco::TrackJet jet2) {
  return jet1.pt() > jet2.pt();
}



//
// constructors and destructor
//
JetCollectionPruner::JetCollectionPruner(const edm::ParameterSet& iConfig)
{
  // define what to produce
  produces<reco::JPTJetCollection>	("jpt"		);
  produces<reco::CaloJetCollection>	("calojet"	);
  produces<reco::PFJetCollection>	("pfjet"	);
  produces<reco::TrackJetCollection>	("trkjet"	);

  // get the input collection of electrons
  inputUncorrectedJPTJetCollection_	= iConfig.getParameter<edm::InputTag>	("inputUncorrectedJPTJetCollection"   	);
  inputUncorrectedPFJetCollection_	= iConfig.getParameter<edm::InputTag>	("inputUncorrectedPFJetCollection" 	);
  inputUncorrectedTrkJetCollection_     = iConfig.getParameter<edm::InputTag>	("inputUncorrectedTrkJetCollection"	); 
  uncorrectedJPTJetPtCut_       	= iConfig.getParameter<double>       	("uncorrectedJPTJetPtCut"            	);
  uncorrectedPFJetPtCut_        	= iConfig.getParameter<double>       	("uncorrectedPFJetPtCut"          	);
  uncorrectedTrkJetPtCut_       	= iConfig.getParameter<double>    	("uncorrectedTrkJetPtCut"		);
}

void JetCollectionPruner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::JPTJetCollection> jetsHandle;
  try {
    iEvent.getByLabel(inputUncorrectedJPTJetCollection_, jetsHandle);
  }
  catch ( cms::Exception& ex ) {
    edm::LogError("JetCollectionPrunerError") << "Error! can't get the input recoJets";
  }

  edm::Handle<reco::PFJetCollection>    pfjetsHandle;
  edm::Handle<reco::TrackJetCollection> trkjetsHandle;
  iEvent.getByLabel(inputUncorrectedPFJetCollection_, pfjetsHandle);
  iEvent.getByLabel(inputUncorrectedTrkJetCollection_, trkjetsHandle);
  
     
  // new collection to put in the event after cleaning
  std::auto_ptr<reco::JPTJetCollection>		outJPTcollection	(new reco::JPTJetCollection);
  std::auto_ptr<reco::CaloJetCollection>	outJETcollection	(new reco::CaloJetCollection);
  std::auto_ptr<reco::PFJetCollection>		outPFJetCollection	(new reco::PFJetCollection);
  std::auto_ptr<reco::TrackJetCollection> 	outTrkJetCollection 	(new reco::TrackJetCollection);

  for (reco::JPTJetCollection::const_iterator jet_it = jetsHandle->begin(); jet_it != jetsHandle->end(); jet_it++)
    {
      if(jet_it->pt() > uncorrectedJPTJetPtCut_)
	outJPTcollection->push_back(*jet_it);
    }

  //sort the jets using the above predicate.
  //they should be sorted, but being paranoid
  std::sort(outJPTcollection->begin(), outJPTcollection->end(), sortJPTJetsByPt);

  for (reco::JPTJetCollection::const_iterator jet_it = outJPTcollection->begin(); jet_it != outJPTcollection->end(); jet_it++)
    {
      const reco::CaloJet& tmp_calojet = dynamic_cast<const reco::CaloJet&>(*(jet_it->getCaloJetRef()));
      outJETcollection->push_back(tmp_calojet);
    }


  for(reco::PFJetCollection::const_iterator jet_it = pfjetsHandle->begin();
      jet_it != pfjetsHandle->end(); jet_it++) {
       
    if(jet_it->pt() > uncorrectedPFJetPtCut_ )
      outPFJetCollection->push_back(*jet_it);
  }

  for(reco::TrackJetCollection::const_iterator jet_it = trkjetsHandle->begin();
      jet_it != trkjetsHandle->end(); jet_it++) {

    if(jet_it->pt() > uncorrectedTrkJetPtCut_ )
      outTrkJetCollection->push_back(*jet_it);

  }
  
  std::sort(outPFJetCollection->begin(),  outPFJetCollection->end(),  sortPFJetsByPt);
  std::sort(outTrkJetCollection->begin(), outTrkJetCollection->end(), sortTrkJetsByPt);


  // put the cleaned collection in the event
  iEvent.put(outJPTcollection,		"jpt"		);
  iEvent.put(outJETcollection,		"calojet"	);
  iEvent.put(outPFJetCollection,	"pfjet"		);
  iEvent.put(outTrkJetCollection,       "trkjet"        );
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

