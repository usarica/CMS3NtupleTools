// -*- C++ -*-
//
// Package:    GenJetMaker
// Class:      GenJetMaker
// 
/**\class GenJetMaker GenJetMaker.cc CMS2/NtupleMaker/src/GenJetMaker.cc

 Description: <one line class summary>
 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Sanjay Padhi
//         Created:  Thu Aug 21 15:47:53 CEST 2008
// $Id: GenJetMaker.cc,v 1.1 2008/08/22 07:49:35 spadhi Exp $
//
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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CMS2/NtupleMaker/interface/GenJetMaker.h"

typedef math::XYZTLorentzVector LorentzVector;

// constructors and destructor
//
GenJetMaker::GenJetMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  produces<unsigned int> ("evtngenjets").setBranchAlias("evt_ngenjets"); // number of jets
  produces<std::vector<LorentzVector> > ("genjetsp4").setBranchAlias("genjets_p4"); // p4 of the jet

  // parameters from configuration
  genJetsInputTag = iConfig.getParameter<edm::InputTag>("genJetsInputTag");
}


GenJetMaker::~GenJetMaker()
{
}

// member functions
// ------------ method called to produce the data  ------------
void
GenJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get GenJet Collection
  edm::Handle<edm::View<reco::GenJet> > genJets;
  iEvent.getByLabel(genJetsInputTag, genJets);
  
  // create containers
  std::auto_ptr<unsigned int> evt_ngenjets(new unsigned int(genJets->size()));
  std::auto_ptr<std::vector<LorentzVector> > vector_genjets_p4(new std::vector<LorentzVector>);
  
  // loop over jets and fill containers
  edm::View<reco::GenJet>::const_iterator jetsEnd = genJets->end();
  for ( edm::View<reco::GenJet>::const_iterator jet = genJets->begin();
        jet != jetsEnd;
        ++jet) {
    vector_genjets_p4->push_back(jet->p4());
  }

  // put containers into event
  iEvent.put(evt_ngenjets, "evtngenjets");
  iEvent.put(vector_genjets_p4, "genjetsp4");
}

void 
GenJetMaker::beginJob(const edm::EventSetup&)
{
}

void 
GenJetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenJetMaker);
