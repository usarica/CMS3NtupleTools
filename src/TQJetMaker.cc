// -*- C++ -*-
//
// Package:    TQJetMaker
// Class:      TQJetMaker
// 
/**\class TQJetMaker TQJetMaker.cc CMS2/NtupleMaker/src/TQJetMaker.cc

Description: copy additional TQAF jet variables in simple data structures into the EDM event tree

 Implementation:
     - take TQAF jets
     - extract and fill variables
*/
//
// Original Author:  pts/4
// Thu Jun 12 22:55:46 UTC 2008
// $Id: TQJetMaker.cc,v 1.1 2008/06/12 23:40:10 gutsche Exp $
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
#include "CMS2/NtupleMaker/interface/TQJetMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;

//
// class decleration
//

//
// constructors and destructor
//
TQJetMaker::TQJetMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  produces<std::vector<int> > ("jetstqgenPartonid").setBranchAlias("jets_tq_genParton_id"); // TQAF gen parton ID
  produces<std::vector<int> > ("jetstqgenPartonMotherid").setBranchAlias("jets_tq_genPartonMother_id"); // TQAF gen parton mother ID
  produces<std::vector<int> > ("jetstqpartonFlavour").setBranchAlias("jets_tq_partonFlavour"); // TQAF parton flavour

  produces<std::vector<float> > ("jetstqnoCorrF").setBranchAlias("jets_tq_noCorrF"); // TQAF corr for corrjet->nocorr
  produces<std::vector<float> > ("jetstqudsCorrF").setBranchAlias("jets_tq_udsCorrF"); // TQAF light quark flavor corr
  produces<std::vector<float> > ("jetstqgluCorrF").setBranchAlias("jets_tq_gluCorrF"); // TQAF gluon  corr
  produces<std::vector<float> > ("jetstqcCorrF").setBranchAlias("jets_tq_cCorrF"); // TQAF c quark flavor corr
  produces<std::vector<float> > ("jetstqbCorrF").setBranchAlias("jets_tq_bCorrF"); // TQAF b quark flavor corr
  produces<std::vector<float> > ("jetstqjetCharge").setBranchAlias("jets_tq_jetCharge"); // TQAF jet charge

  produces<std::vector<LorentzVector> > ("jetstqgenPartonp4").setBranchAlias("jets_tq_genParton_p4"); // TQAF gen parton p4
  produces<std::vector<LorentzVector> > ("jetstqgenPartonMotherp4").setBranchAlias("jets_tq_genPartonMother_p4"); // TQAF gen parton mother p4
  produces<std::vector<LorentzVector> > ("jetstqgenJetp4").setBranchAlias("jets_tq_genJet_p4"); // TQAF gen jet p4
  produces<std::vector<LorentzVector> > ("jetstqjetp4").setBranchAlias("jets_tq_jet_p4"); // TQAF jet p4

  // parameters from configuration
  tqJetsInputTag = iConfig.getParameter<edm::InputTag>("tqJetsInputTag");

}


TQJetMaker::~TQJetMaker()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TQJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get jet collection
  edm::Handle<std::vector<TopJet> > tqJetsHandle;
  iEvent.getByLabel(tqJetsInputTag, tqJetsHandle);

  // create containers
  std::auto_ptr<std::vector<int> > vector_jets_tq_genParton_id(new std::vector<int>);
  std::auto_ptr<std::vector<int> > vector_jets_tq_genPartonMother_id(new std::vector<int>);
  std::auto_ptr<std::vector<int> > vector_jets_tq_partonFlavour(new std::vector<int>);

  std::auto_ptr<std::vector<float> > vector_jets_tq_noCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_tq_udsCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_tq_gluCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_tq_cCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_tq_bCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_tq_jetCharge(new std::vector<float>);

  std::auto_ptr<std::vector<LorentzVector> > vector_jets_tq_genParton_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_jets_tq_genPartonMother_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_jets_tq_genJet_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_jets_tq_jet_p4(new std::vector<LorentzVector>);

  // loop over jets and fill containers
  std::vector<TopJet>::const_iterator tqJetsEnd = tqJetsHandle->end(); 
  for ( std::vector<TopJet>::const_iterator tqJet = tqJetsHandle->begin();
	tqJet != tqJetsEnd; 
	++tqJet) {

    const reco::GenParticle *mother = MCUtilities::motherID(tqJet->getGenParton());

    vector_jets_tq_genParton_id->push_back(tqJet->getGenParton().pdgId());
    vector_jets_tq_genPartonMother_id->push_back(mother ? mother->pdgId() : 0);
    vector_jets_tq_partonFlavour->push_back(tqJet->getPartonFlavour());

    vector_jets_tq_noCorrF->push_back(tqJet->getNoCorrF());
    vector_jets_tq_udsCorrF->push_back(tqJet->getUdsCorrF());
    vector_jets_tq_gluCorrF->push_back(tqJet->getGluCorrF());
    vector_jets_tq_cCorrF->push_back(tqJet->getCCorrF());
    vector_jets_tq_bCorrF->push_back(tqJet->getBCorrF());
    vector_jets_tq_jetCharge->push_back(tqJet->getJetCharge());

    vector_jets_tq_genParton_p4->push_back(tqJet->getGenParton().p4());
    vector_jets_tq_genPartonMother_p4->push_back(mother ? mother->p4() : LorentzVector(0,0,0,0) );
    vector_jets_tq_genJet_p4->push_back(tqJet->getGenJet().p4());
    vector_jets_tq_jet_p4->push_back(tqJet->p4());
  }

  // put containers into event
  iEvent.put(vector_jets_tq_genParton_id, "jetstqgenPartonid");
  iEvent.put(vector_jets_tq_genPartonMother_id, "jetstqgenPartonMotherid");
  iEvent.put(vector_jets_tq_partonFlavour, "jetstqpartonFlavour");

  iEvent.put(vector_jets_tq_noCorrF, "jetstqnoCorrF");
  iEvent.put(vector_jets_tq_udsCorrF, "jetstqudsCorrF");
  iEvent.put(vector_jets_tq_gluCorrF, "jetstqgluCorrF");
  iEvent.put(vector_jets_tq_cCorrF, "jetstqcCorrF");
  iEvent.put(vector_jets_tq_bCorrF, "jetstqbCorrF");
  iEvent.put(vector_jets_tq_jetCharge, "jetstqjetCharge");

  iEvent.put(vector_jets_tq_genParton_p4, "jetstqgenPartonp4");
  iEvent.put(vector_jets_tq_genPartonMother_p4, "jetstqgenPartonMotherp4");
  iEvent.put(vector_jets_tq_genJet_p4, "jetstqgenJetp4");
  iEvent.put(vector_jets_tq_jet_p4, "jetstqjetp4");
}

// ------------ method called once each job just before starting event loop  ------------
void 
TQJetMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TQJetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TQJetMaker);
