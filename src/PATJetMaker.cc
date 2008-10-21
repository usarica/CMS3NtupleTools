// -*- C++ -*-
//
// Package:    PATJetMaker
// Class:      PATJetMaker
// 
/**\class PATJetMaker PATJetMaker.cc CMS2/NtupleMaker/src/PATJetMaker.cc

Description: copy additional PAT jet variables in simple data structures into the EDM event tree

 Implementation:
     - take PAT jets
     - extract and fill variables
*/
//
// Original Author:  pts/4
// Thu Jun 12 22:55:46 UTC 2008
// $Id: PATJetMaker.cc,v 1.1 2008/10/21 07:26:20 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/PATJetMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Flags.h"

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
PATJetMaker::PATJetMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  produces<std::vector<int> > ("jetspatgenPartonid").setBranchAlias("jets_pat_genParton_id"); // PAT gen parton ID
  produces<std::vector<int> > ("jetspatgenPartonMotherid").setBranchAlias("jets_pat_genPartonMother_id"); // PAT gen parton mother ID
  produces<std::vector<int> > ("jetspatpartonFlavour").setBranchAlias("jets_pat_partonFlavour"); // PAT parton flavour
  produces<std::vector<uint32_t> >("jetspatflag").setBranchAlias("jets_pat_flag"); //PAT flag

  produces<std::vector<float> > ("jetspatnoCorrF").setBranchAlias("jets_pat_noCorrF"); // PAT corr for corrjet->nocorr
  produces<std::vector<float> > ("jetspatudsCorrF").setBranchAlias("jets_pat_udsCorrF"); // PAT light quark flavor corr
  produces<std::vector<float> > ("jetspatgluCorrF").setBranchAlias("jets_pat_gluCorrF"); // PAT gluon  corr
  produces<std::vector<float> > ("jetspatcCorrF").setBranchAlias("jets_pat_cCorrF"); // PAT c quark flavor corr
  produces<std::vector<float> > ("jetspatbCorrF").setBranchAlias("jets_pat_bCorrF"); // PAT b quark flavor corr
  produces<std::vector<float> > ("jetspatjetCharge").setBranchAlias("jets_pat_jetCharge"); // PAT jet charge

  produces<std::vector<LorentzVector> > ("jetspatgenPartonp4").setBranchAlias("jets_pat_genParton_p4"); // PAT gen parton p4
  produces<std::vector<LorentzVector> > ("jetspatgenPartonMotherp4").setBranchAlias("jets_pat_genPartonMother_p4"); // PAT gen parton mother p4
  produces<std::vector<LorentzVector> > ("jetspatgenJetp4").setBranchAlias("jets_pat_genJet_p4"); // PAT gen jet p4
  produces<std::vector<LorentzVector> > ("jetspatjetp4").setBranchAlias("jets_pat_jet_p4"); // PAT jet p4

  // parameters from configuration
  patJetsInputTag = iConfig.getParameter<edm::InputTag>("patJetsInputTag");

}


PATJetMaker::~PATJetMaker()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PATJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get jet collection
  edm::Handle<std::vector<pat::Jet> > patJetsHandle;
  iEvent.getByLabel(patJetsInputTag, patJetsHandle);

  // create containers
  std::auto_ptr<std::vector<int> > vector_jets_pat_genParton_id(new std::vector<int>);
  std::auto_ptr<std::vector<int> > vector_jets_pat_genPartonMother_id(new std::vector<int>);
  std::auto_ptr<std::vector<int> > vector_jets_pat_partonFlavour(new std::vector<int>);
  std::auto_ptr<std::vector<uint32_t> > vector_jets_pat_flag(new std::vector<uint32_t>);

  std::auto_ptr<std::vector<float> > vector_jets_pat_noCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_pat_udsCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_pat_gluCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_pat_cCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_pat_bCorrF(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jets_pat_jetCharge(new std::vector<float>);

  std::auto_ptr<std::vector<LorentzVector> > vector_jets_pat_genParton_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_jets_pat_genPartonMother_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_jets_pat_genJet_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<LorentzVector> > vector_jets_pat_jet_p4(new std::vector<LorentzVector>);

  // loop over jets and fill containers
  std::vector<pat::Jet>::const_iterator patJetsEnd = patJetsHandle->end(); 
  for ( std::vector<pat::Jet>::const_iterator patJet = patJetsHandle->begin();
	patJet != patJetsEnd; 
	++patJet) {

    reco::GenParticle genParton = patJet->genParton() ? *patJet->genParton() : 
      reco::GenParticle(0, reco::Particle::LorentzVector(0, 0, 0, 0), reco::Particle::Point(0,0,0), 0, 0, true);

    const reco::GenParticle *mother = MCUtilities::motherID(genParton);

    vector_jets_pat_genParton_id->push_back(genParton.pdgId());
    vector_jets_pat_genPartonMother_id->push_back(mother ? mother->pdgId() : 0);
    vector_jets_pat_partonFlavour->push_back(patJet->partonFlavour());
    vector_jets_pat_flag->push_back(patJet->status());

    vector_jets_pat_noCorrF->push_back(patJet->correctionFactor(pat::Jet::NoCorrection));
    vector_jets_pat_udsCorrF->push_back(patJet->correctionFactor(pat::Jet::udsCorrection));
    vector_jets_pat_gluCorrF->push_back(patJet->correctionFactor(pat::Jet::gCorrection));
    vector_jets_pat_cCorrF->push_back(patJet->correctionFactor(pat::Jet::cCorrection));
    vector_jets_pat_bCorrF->push_back(patJet->correctionFactor(pat::Jet::bCorrection));
    vector_jets_pat_jetCharge->push_back(patJet->jetCharge());

    vector_jets_pat_genParton_p4->push_back(genParton.p4());
    vector_jets_pat_genPartonMother_p4->push_back(mother ? mother->p4() : LorentzVector(0,0,0,0) );
    LorentzVector genJetP4 = patJet->genJet() ? patJet->genJet()->p4() : LorentzVector(0, 0, 0, 0);
    vector_jets_pat_genJet_p4->push_back(genJetP4);
    vector_jets_pat_jet_p4->push_back(patJet->p4());
  }

  // put containers into event
  iEvent.put(vector_jets_pat_genParton_id, "jetspatgenPartonid");
  iEvent.put(vector_jets_pat_genPartonMother_id, "jetspatgenPartonMotherid");
  iEvent.put(vector_jets_pat_partonFlavour, "jetspatpartonFlavour");
  iEvent.put(vector_jets_pat_flag, "jetspatflag");     
  iEvent.put(vector_jets_pat_noCorrF, "jetspatnoCorrF");
  iEvent.put(vector_jets_pat_udsCorrF, "jetspatudsCorrF");
  iEvent.put(vector_jets_pat_gluCorrF, "jetspatgluCorrF");
  iEvent.put(vector_jets_pat_cCorrF, "jetspatcCorrF");
  iEvent.put(vector_jets_pat_bCorrF, "jetspatbCorrF");
  iEvent.put(vector_jets_pat_jetCharge, "jetspatjetCharge");

  iEvent.put(vector_jets_pat_genParton_p4, "jetspatgenPartonp4");
  iEvent.put(vector_jets_pat_genPartonMother_p4, "jetspatgenPartonMotherp4");
  iEvent.put(vector_jets_pat_genJet_p4, "jetspatgenJetp4");
  iEvent.put(vector_jets_pat_jet_p4, "jetspatjetp4");
}

// ------------ method called once each job just before starting event loop  ------------
void 
PATJetMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATJetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATJetMaker);
