// -*- C++ -*-
//
// Package:    JetMaker
// Class:      JetMaker
// 
/**\class JetMaker JetMaker.cc CMS2/NtupleMaker/src/JetMaker.cc

Description: copy reco::CaloJet variables in simple data structures into the EDM event tree

 Implementation:
     - take TQAF jets
     - extract and fill variables
*/
//
// Original Author:  Oliver Gutsche
// Created:  Tue Jun  9 11:07:38 CDT 2008
// $Id: JetMaker.cc,v 1.2 2008/06/12 23:40:37 gutsche Exp $
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
#include "CMS2/NtupleMaker/interface/JetMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;

//
// class decleration
//

//
// constructors and destructor
//
JetMaker::JetMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  produces<unsigned int> ("njets").setBranchAlias("njets"); // number of jets
  produces<std::vector<LorentzVector> >	("jetsp4").setBranchAlias("jets_p4"); // p4 of the jet
  produces<std::vector<LorentzVector> >	("jetsmcp4").setBranchAlias("jets_mc_p4"); // p4 of the matched GenJet
  produces<std::vector<LorentzVector> >	("jetsmcgpp4").setBranchAlias("jets_mc_gp_p4"); // p4 of the matched MC particle
  produces<std::vector<float> >	        ("jetsemFrac").setBranchAlias("jets_emFrac"); // electromagnetic energy fraction
  produces<std::vector<float> >	        ("jetschFrac").setBranchAlias("jets_chFrac"); // charged track energy fraction
  produces<std::vector<float> >	        ("jetsmcemEnergy").setBranchAlias("jets_mc_emEnergy"); // energy of electromagnetic particles of the matched GenJet
  produces<std::vector<float> >	        ("jetsmchadEnergy").setBranchAlias("jets_mc_hadEnergy"); // energy of hadronic particles of the matched GenJet
  produces<std::vector<float> >	        ("jetsmcinvEnergy").setBranchAlias("jets_mc_invEnergy"); // invisible energy of the matched GenJet
  produces<std::vector<float> >	        ("jetsmcotherEnergy").setBranchAlias("jets_mc_otherEnergy"); // other energy (undecayed Sigmas etc.) of the matched GenJet
  produces<std::vector<float> >	        ("jetscor").setBranchAlias("jets_cor"); // energy scale correction
  produces<std::vector<float> >	        ("jetsEMFcor").setBranchAlias("jets_EMFcor"); // energy scale corrections including electromagnetic fraction of jet
  produces<std::vector<int> >           ("jetsmcid").setBranchAlias("jets_mc_id"); // PDG id of the matched MC particle

  // parameters from configuration
  jetsInputTag = iConfig.getParameter<edm::InputTag>("jetsInputTag");
  genJetsInputTag = iConfig.getParameter<edm::InputTag>("genJetsInputTag");
  genParticlesInputTag = iConfig.getParameter<edm::InputTag>("genParticlesInputTag");

}


JetMaker::~JetMaker()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
JetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get jet collection
  edm::Handle<edm::View<reco::CaloJet> > jetsHandle;
  iEvent.getByLabel(jetsInputTag, jetsHandle);

  // get generator jet collection
  edm::Handle<reco::GenJetCollection> genJetsHandle;
  iEvent.getByLabel(genJetsInputTag, genJetsHandle);

  // get MC particle collection
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  iEvent.getByLabel(genParticlesInputTag, genParticlesHandle);

  // create containers
  std::auto_ptr<unsigned int> njets(new unsigned int(jetsHandle->size()));
  std::auto_ptr<std::vector<LorentzVector> > vector_jets_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_jets_mc_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_jets_mc_gp_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<float> > vector_jets_emFrac(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_chFrac(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_mc_emEnergy(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_mc_hadEnergy(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_mc_invEnergy(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_mc_otherEnergy(new std::vector<float>);
//   std::auto_ptr<std::vector<float> > vector_jets_cor(new std::vector<float>);
//   std::auto_ptr<std::vector<float> > vector_jets_EMFcor(new std::vector<float>);
  std::auto_ptr<std::vector<int> > vector_jets_mc_id(new std::vector<int>);

  // loop over jets and fill containers
  edm::View<reco::CaloJet>::const_iterator jetsEnd = jetsHandle->end(); 
  for ( edm::View<reco::CaloJet>::const_iterator jet = jetsHandle->begin();
	jet != jetsEnd; 
	++jet) {
    vector_jets_p4->push_back(jet->p4());
    vector_jets_emFrac->push_back(jet->emEnergyFraction());
    vector_jets_chFrac->push_back(-999.);

    const reco::GenJet* matchedGenJet = MatchUtilities::matchCandToGenJet(*jet,genJetsHandle.product());
    if ( matchedGenJet != 0 ) {
      vector_jets_mc_p4->push_back(matchedGenJet->p4());
      vector_jets_mc_emEnergy->push_back(matchedGenJet->emEnergy());
      vector_jets_mc_hadEnergy->push_back(matchedGenJet->hadEnergy());
      vector_jets_mc_invEnergy->push_back(matchedGenJet->invisibleEnergy());
      vector_jets_mc_otherEnergy->push_back(matchedGenJet->auxiliaryEnergy());
    } else {
      edm::LogWarning("JetMaker") << "Jet " << (jet - jetsHandle->begin()) << " could not be matched to a GetJet!";
      vector_jets_mc_p4->push_back(LorentzVector(0,0,0,0));
      vector_jets_mc_emEnergy->push_back(-999.);
      vector_jets_mc_hadEnergy->push_back(-999.);
      vector_jets_mc_invEnergy->push_back(-999.);
      vector_jets_mc_otherEnergy->push_back(-999.);
    }

    const reco::GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*jet,genParticlesHandle.product());
    if ( matchedGenParticle != 0 ) {
      vector_jets_mc_gp_p4->push_back(matchedGenParticle->p4());
      vector_jets_mc_id->push_back(matchedGenParticle->pdgId());
    } else {
      edm::LogWarning("JetMaker") << "Jet " << (jet - jetsHandle->begin()) << " could not be matched to a GetParticle!";
      vector_jets_mc_gp_p4->push_back(LorentzVector(0,0,0,0));
      vector_jets_mc_id->push_back(0);
    }

//     vector_jets_cor->push_back(-999.);
//     vector_jets_EMFcor->push_back(-999.);

  }

  // put containers into event
  iEvent.put(njets, "njets");
  iEvent.put(vector_jets_p4, "jetsp4");
  iEvent.put(vector_jets_mc_p4,"jetsmcp4");
  iEvent.put(vector_jets_mc_gp_p4,"jetsmcgpp4");
  iEvent.put(vector_jets_emFrac,"jetsemFrac");
  iEvent.put(vector_jets_chFrac,"jetschFrac");
  iEvent.put(vector_jets_mc_emEnergy,"jetsmcemEnergy");
  iEvent.put(vector_jets_mc_hadEnergy,"jetsmchadEnergy");
  iEvent.put(vector_jets_mc_invEnergy,"jetsmcinvEnergy");
  iEvent.put(vector_jets_mc_otherEnergy,"jetsmcotherEnergy");
//   iEvent.put(vector_jets_cor,"jetscor");
//   iEvent.put(vector_jets_EMFcor,"jetsEMFcor");
  iEvent.put(vector_jets_mc_id,"jetsmcid");
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetMaker);
