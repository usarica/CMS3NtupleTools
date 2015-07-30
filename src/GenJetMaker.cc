// -*- C++ -*-
//
// Package:    GenJetMaker
// Class:      GenJetMaker
// 
/**\class GenJetMaker GenJetMaker.cc CMS3/NtupleMaker/src/GenJetMaker.cc

   Description: <one line class summary>
   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Sanjay Padhi
//         Created:  Thu Aug 21 15:47:53 CEST 2008
// $Id: GenJetMaker.cc,v 1.7 2013/02/18 17:23:22 dalfonso Exp $
//
//

// system include files
#include <memory>
#include <vector>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "CMS3/NtupleMaker/interface/GenJetMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;

bool sortByPt(LorentzVector jet1, LorentzVector jet2) {
  return jet1.pt() > jet2.pt();
}

// constructors and destructor
//
GenJetMaker::GenJetMaker(const edm::ParameterSet& iConfig)
{

  aliaspostfix_ = iConfig.getUntrackedParameter<std::string>("aliasPostfix");

  // product of this EDProducer
  produces<unsigned int>                ("evtngenjets"+aliaspostfix_).setBranchAlias("evt_ngenjets"+aliaspostfix_); // number of jets
  produces<std::vector<LorentzVector> > ("genjetsp4"+aliaspostfix_).setBranchAlias("genjets_p4"+aliaspostfix_); // p4 of the jet

  // parameters from configuration
  genJetsInputTag = iConfig.getParameter<edm::InputTag>("genJetsInputTag");
  genJetMinPtCut  = iConfig.getParameter<double>       ("genJetMinPtCut" );
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

  if ( !genJets.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve gen jets collection";
    edm::LogInfo("OutputInfo") << " GenJetMaker cannot continue...!";
    std::cout << "GenJetMaker cannot continue" << std::endl;
    throw cms::Exception("GenJetMaker cannot continue!  Could not load gen jets collection.");
  }
  
  // create containers
  std::auto_ptr<unsigned int>                evt_ngenjets     (new unsigned int(genJets->size()) );
  std::auto_ptr<std::vector<LorentzVector> > vector_genjets_p4(new std::vector<LorentzVector>    );
  
  // loop over jets and fill containers
  edm::View<reco::GenJet>::const_iterator jetsEnd = genJets->end();
  for ( edm::View<reco::GenJet>::const_iterator jet = genJets->begin(); jet != jetsEnd; ++jet) {

    if( jet->pt() < genJetMinPtCut ) continue;

    vector_genjets_p4->push_back( LorentzVector( jet->p4() ) );
  }

  // sort gen jets by pt
  std::sort(vector_genjets_p4->begin(), vector_genjets_p4->end(), sortByPt);

  // put containers into event
  iEvent.put(evt_ngenjets     , "evtngenjets"+aliaspostfix_);
  iEvent.put(vector_genjets_p4, "genjetsp4"+aliaspostfix_);

}

void GenJetMaker::beginJob()
{
}

void GenJetMaker::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenJetMaker);
