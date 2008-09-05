//-*- C++ -*-
//
// Package:    AllTrkJetMaker
// Class:      AllTrkJetMaker.cc
//
/**\class AllTrkJetMaker AllTrkJetMaker.cc CMS2/NtupleMaker/src/AllTrkJetMaker.cc

Description: Dumps the TrkJet contents into the ntuple
Implementation:
*/
//
// Original Author:  Sanjay Padhi
//         Created:  Mon Jun 23 03:57:47 CEST 2008
// $Id: AllTrkJetMaker.cc,v 1.1 2008/09/05 00:22:16 fkw Exp $
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/AllTrkJetMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;


//
// constructors and destructor
//
AllTrkJetMaker::AllTrkJetMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  produces<unsigned int> ("evtntrkjets").setBranchAlias("evt_nalltrkjets"); // number of jets
  produces<std::vector<LorentzVector> > ("trkjetsp4").setBranchAlias("alltrkjets_p4"); // p4 of the jet

  // parameters from configuration
  trkJetsInputTag = iConfig.getParameter<edm::InputTag>("trkJetsInputTag");

}


AllTrkJetMaker::~AllTrkJetMaker()
{

}

void
AllTrkJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get TrkJet collection
  edm::Handle<edm::View<reco::BasicJet> > trkJets;
  iEvent.getByLabel(trkJetsInputTag, trkJets);
  
  // create containers
  std::auto_ptr<unsigned int> evt_ntrkjets(new unsigned int(trkJets->size()));
  std::auto_ptr<std::vector<LorentzVector> > vector_trkjets_p4(new std::vector<LorentzVector>);
  
  // loop over jets and fill containers
  edm::View<reco::BasicJet>::const_iterator jetsEnd = trkJets->end();
  for ( edm::View<reco::BasicJet>::const_iterator jet = trkJets->begin();
        jet != jetsEnd;
        ++jet) {
    vector_trkjets_p4->push_back(jet->p4());
  }

  // put containers into event
  iEvent.put(evt_ntrkjets, "evtntrkjets");
  iEvent.put(vector_trkjets_p4, "trkjetsp4");

}


void
AllTrkJetMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
AllTrkJetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(AllTrkJetMaker);
