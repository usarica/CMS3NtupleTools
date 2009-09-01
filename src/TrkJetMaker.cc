//-*- C++ -*-
//
// Package:    TrkJetMaker
// Class:      TrkJetMaker.cc
//
/**\class TrkJetMaker TrkJetMaker.cc CMS2/NtupleMaker/src/TrkJetMaker.cc

Description: Dumps the TrkJet contents into the ntuple
Implementation:
*/
//
// Original Author:  Sanjay Padhi
//         Created:  Mon Jun 23 03:57:47 CEST 2008
// $Id: TrkJetMaker.cc,v 1.2 2009/09/01 11:49:00 fgolf Exp $
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
#include "CMS2/NtupleMaker/interface/TrkJetMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;

bool sortJetsByPt(LorentzVector jet1, LorentzVector jet2) {
  return jet1.pt() > jet2.pt();
}

//
// constructors and destructor
//
TrkJetMaker::TrkJetMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  produces<unsigned int> ("evtntrkjets").setBranchAlias("evt_ntrkjets"); // number of jets
  produces<std::vector<LorentzVector> > ("trkjetsp4").setBranchAlias("trkjets_p4"); // p4 of the jet

  // parameters from configuration
  trkJetsInputTag = iConfig.getParameter<edm::InputTag>("trkJetsInputTag");
  trkJetPtCut_    = iConfig.getParameter<double>       ("trkJetPtCut"    );
}


TrkJetMaker::~TrkJetMaker()
{

}

void
TrkJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get TrkJet collection
  edm::Handle<edm::View<reco::BasicJet> > trkJets;
  iEvent.getByLabel(trkJetsInputTag, trkJets);
  
  // create containers
  std::auto_ptr<unsigned int> evt_ntrkjets(new unsigned int(trkJets->size()));
  std::auto_ptr<std::vector<LorentzVector> > vector_trkjets_p4(new std::vector<LorentzVector>);
  
  // loop over jets and fill containers
  edm::View<reco::BasicJet>::const_iterator jetsEnd = trkJets->end();
  for ( edm::View<reco::BasicJet>::const_iterator jet = trkJets->begin(); jet != jetsEnd; ++jet) {

    if( jet->p4().pt() < trkJetPtCut_ ) continue;

    vector_trkjets_p4->push_back(jet->p4());
  }

  std::sort( vector_trkjets_p4->begin(), vector_trkjets_p4->end(), sortJetsByPt );

  // put containers into event
  iEvent.put(evt_ntrkjets, "evtntrkjets");
  iEvent.put(vector_trkjets_p4, "trkjetsp4");

}


void
TrkJetMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TrkJetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrkJetMaker);
