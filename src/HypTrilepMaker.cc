// -*- C++ -*-
//
// Package:    HypTrilepMaker
// Class:      HypTrilepMaker
// 
/**\class HypTrilepMaker HypTrilepMaker.cc CMS2/NtupleMaker/src/HypTrilepMaker.cc

Description: create trilepton hypothesis branches

Implementation:
- combine muons and electrons after preselection
- correct jets and store index vectors
- correct met

Hypothesis:
- lepton type: 1:m+, -1:m-, 2:e+, -2:e-
- trilepton 20 buckets:
e+e+e+: 0
e+e+e-: 1
e+e-e-: 2
e-e-e-: 3
m+m+m+: 4
m+m+m-: 5
m+m-m-: 6
m-m-m-: 7
m+e+e+: 8
m-e+e+: 9
m-e+e-: 10
m+e+e-: 11
m-e-e-: 12
m+e-e-: 13
m+m+e+: 14
m+m+e-: 15
m+m-e-: 16
m+m-e+: 17
m-m-e-: 18
m-m-e+: 19
- trilepton hypothesis first/second/third follows bucket definition. In the case of the same type, ordering follows descending index number within muon and electron collection.

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: HypTrilepMaker.cc,v 1.2 2008/06/25 00:36:12 gutsche Exp $
//
//

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/HypTrilepMaker.h"


#include "CMS2/NtupleMaker/interface/MatchUtilities.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//
// class decleration
//

//
// constructors and destructor
//
HypTrilepMaker::HypTrilepMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  // 
  // trilepton hyptothesis
  //
  produces<std::vector<unsigned int> > ("hyptrilepbucket").setBranchAlias("hyp_trilep_bucket");            // trilepton bucket
  produces<std::vector<int> >          ("hyptrilepfirsttype").setBranchAlias("hyp_trilep_first_type");     // type of the first lepton in the trilepton hypothesis (1: muon, 2: electron)
  produces<std::vector<unsigned int> > ("hyptrilepfirstindex").setBranchAlias("hyp_trilep_first_index");   // index of first lepton in lepton collection
  produces<std::vector<int> >          ("hyptrilepsecondtype").setBranchAlias("hyp_trilep_second_type");   // type of the second lepton in the trilepton hypothesis (1: muon, 2: electron)
  produces<std::vector<unsigned int> > ("hyptrilepsecondindex").setBranchAlias("hyp_trilep_second_index"); // index of second lepton in lepton collection
  produces<std::vector<int> >          ("hyptrilepthirdtype").setBranchAlias("hyp_trilep_third_type");     // type of the third lepton in the trilepton hypothesis (1: muon, 2: electron)
  produces<std::vector<unsigned int> > ("hyptrilepthirdindex").setBranchAlias("hyp_trilep_third_index");   // index of third lepton in lepton collection

  // parameters from configuration
  muonsInputTag = iConfig.getParameter<edm::InputTag>("muonsInputTag");
  electronsInputTag = iConfig.getParameter<edm::InputTag>("electronsInputTag");

}


HypTrilepMaker::~HypTrilepMaker()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HypTrilepMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //input collections

  // muon charge
  edm::InputTag mus_charge_tag(muonsInputTag.label(),"muscharge");
  edm::Handle<std::vector<int> > mus_charge_h;
  iEvent.getByLabel(mus_charge_tag, mus_charge_h);
  const std::vector<int> *mus_charge = mus_charge_h.product();

  // electron charge
  edm::InputTag els_charge_tag(electronsInputTag.label(),"elscharge");
  edm::Handle<std::vector<int> > els_charge_h;
  iEvent.getByLabel(els_charge_tag, els_charge_h);
  const std::vector<int> *els_charge = els_charge_h.product();

  // output collections
  std::auto_ptr<std::vector<unsigned int> > vector_hyp_trilep_bucket(new std::vector<unsigned int>);
  std::auto_ptr<std::vector<int> > vector_hyp_trilep_first_type(new std::vector<int>);
  std::auto_ptr<std::vector<unsigned int> > vector_hyp_trilep_first_index(new std::vector<unsigned int>);
  std::auto_ptr<std::vector<int> > vector_hyp_trilep_second_type(new std::vector<int>);
  std::auto_ptr<std::vector<unsigned int> > vector_hyp_trilep_second_index(new std::vector<unsigned int>);
  std::auto_ptr<std::vector<int> > vector_hyp_trilep_third_type(new std::vector<int>);
  std::auto_ptr<std::vector<unsigned int> > vector_hyp_trilep_third_index(new std::vector<unsigned int>);

  // processed trilepton candidates 
  // ordered: first by type (1st muon, 2nd electron), then by index number
  std::vector<unsigned int> cand;

  // array for sorts
  unsigned int sorter[3] = {0,0,0};


  // number of electrons
  unsigned int evt_nels = els_charge->size();

  // number of electrons
  unsigned int evt_nmus = mus_charge->size();

  // m
  for (unsigned int firstMuon = 0; firstMuon < evt_nmus; ++firstMuon) {
    // m
    for (unsigned int secondMuon = 0; secondMuon < evt_nmus; ++secondMuon) {
      if ( secondMuon == firstMuon ) continue;
      // m
      for (unsigned int thirdMuon = 0; thirdMuon < evt_nmus; ++thirdMuon) {
	if ( thirdMuon == firstMuon || thirdMuon == secondMuon ) continue;
	sorter[0] = firstMuon;
	sorter[1] = secondMuon;
	sorter[2] = thirdMuon;
	std::sort(sorter,sorter+3);
	unsigned int candIndex = encodeTrileptonCandidate(sorter[0],
							  sorter[1],
							  sorter[2]);
	if ( find(cand.begin(),cand.end(), candIndex) == cand.end() ) {
	  cand.push_back(candIndex);
	  int charge = mus_charge->at(sorter[0]) + mus_charge->at(sorter[1]) + mus_charge->at(sorter[2]);
	  if ( charge == -3 ) {
	    // m-m-m-
	    vector_hyp_trilep_bucket->push_back(7);
	  } else if ( charge == -1) {
	    // m+m-m-
	    vector_hyp_trilep_bucket->push_back(6);
	  } else if ( charge == 1 ) {
	    // m+m+m-
	    vector_hyp_trilep_bucket->push_back(5);
	  } else if ( charge == 3 ) {
	    // m+m+m+
	    vector_hyp_trilep_bucket->push_back(4);
	  } else {
	    edm::LogError("HypTrilepMaker") << "combineTriLeptons mmm : charge combination could not be identified!!!";
	  }
	  vector_hyp_trilep_first_type->push_back(mus_charge->at(sorter[0]));
	  vector_hyp_trilep_first_index->push_back(sorter[0]);
	  vector_hyp_trilep_second_type->push_back(mus_charge->at(sorter[1]));
	  vector_hyp_trilep_second_index->push_back(sorter[1]);
	  vector_hyp_trilep_third_type->push_back(mus_charge->at(sorter[2]));
	  vector_hyp_trilep_third_index->push_back(sorter[2]);
	}
      }
      // e
      for (unsigned int thirdElectron = 0; thirdElectron < evt_nels; ++thirdElectron) {
	// order muon indices
	sorter[0] = firstMuon;
	sorter[1] = secondMuon;
	std::sort(sorter,sorter+2);
	// add electron
	sorter[2] = thirdElectron;
	unsigned int candIndex = encodeTrileptonCandidate(sorter[0],
							  sorter[1],
							  sorter[2]);
	if ( find(cand.begin(),cand.end(), candIndex) == cand.end() ) {
	  cand.push_back(candIndex);
	  int charge = mus_charge->at(sorter[0]) + mus_charge->at(sorter[1]);
	  if ( charge == -2 ) {
	    if ( els_charge->at(sorter[2]) < 0 ) {
	      // m-m-e-
	      vector_hyp_trilep_bucket->push_back(18);
	    } else {
	      // m-m-e+
	      vector_hyp_trilep_bucket->push_back(19);
	    }
	  } else if ( charge == 0) {
	    if ( els_charge->at(sorter[2]) < 0 ) {
	      // m+m-e-
	      vector_hyp_trilep_bucket->push_back(16);
	    } else {
	      // m+m-e+
	      vector_hyp_trilep_bucket->push_back(17);
	    }
	  } else if ( charge == 2 ) {
	    if ( els_charge->at(sorter[2]) < 0 ) {
	      // m+m+e-
	      vector_hyp_trilep_bucket->push_back(15);
	    } else {
	      // m+m+e+
	      vector_hyp_trilep_bucket->push_back(14);
	    }
	  } else {
	    edm::LogError("HypTrilepMaker") << "combineTriLeptons mme: charge combination could not be identified!!!";
	  }
	  vector_hyp_trilep_first_type->push_back(mus_charge->at(sorter[0]));
	  vector_hyp_trilep_first_index->push_back(sorter[0]);
	  vector_hyp_trilep_second_type->push_back(mus_charge->at(sorter[1]));
	  vector_hyp_trilep_second_index->push_back(sorter[1]);
	  vector_hyp_trilep_third_type->push_back(els_charge->at(sorter[2])*2);
	  vector_hyp_trilep_third_index->push_back(sorter[2]);
	}
      }
    }
    // e
    for (unsigned int secondElectron = 0; secondElectron < evt_nels; ++secondElectron) {
      // e
      for (unsigned int thirdElectron = 0; thirdElectron < evt_nels; ++thirdElectron) {
	if ( thirdElectron == secondElectron ) continue;
	// order electron indices
	sorter[0] = secondElectron;
	sorter[1] = thirdElectron;
	std::sort(sorter,sorter+2);
	// add muon in first place
	sorter[2] = sorter[1];
	sorter[1] = sorter[0];
	sorter[0] = firstMuon;
	unsigned int candIndex = encodeTrileptonCandidate(sorter[0],
							  sorter[1],
							  sorter[2]);
	if ( find(cand.begin(),cand.end(), candIndex) == cand.end() ) {
	  cand.push_back(candIndex);
	  // check for charge to identify bucket number
	  int charge = els_charge->at(sorter[1]) + els_charge->at(sorter[2]);
	  if ( charge == -2 ) {
	    if ( mus_charge->at(sorter[0]) < 0 ) {
	      // e-e-m-
	      vector_hyp_trilep_bucket->push_back(12);
	    } else {
	      // e-e-m+
	      vector_hyp_trilep_bucket->push_back(13);
	    }
	  } else if ( charge == 0) {
	    if ( mus_charge->at(sorter[0]) < 0 ) {
	      // e+e-m-
	      vector_hyp_trilep_bucket->push_back(10);
	    } else {
	      // e+e-m+
	      vector_hyp_trilep_bucket->push_back(11);
	    }
	  } else if ( charge == 2 ) {
	    if ( mus_charge->at(sorter[0]) < 0 ) {
	      // e+e+m-
	      vector_hyp_trilep_bucket->push_back(9);
	    } else {
	      // e+e+m+
	      vector_hyp_trilep_bucket->push_back(8);
	    }
	  } else {
	    edm::LogError("HypTrilepMaker") << "combineTriLeptons mee: charge combination could not be identified!!!";
	  }
	  vector_hyp_trilep_first_type->push_back(mus_charge->at(sorter[0]));
	  vector_hyp_trilep_first_index->push_back(sorter[0]);
	  vector_hyp_trilep_second_type->push_back(els_charge->at(sorter[1])*2);
	  vector_hyp_trilep_second_index->push_back(sorter[1]);
	  vector_hyp_trilep_third_type->push_back(els_charge->at(sorter[2])*2);
	  vector_hyp_trilep_third_index->push_back(sorter[2]);
          

	}
      }
    }
  }
  
  // e
  for (unsigned int firstElectron = 0; firstElectron < evt_nels; ++firstElectron) {
    // e
    for (unsigned int secondElectron = 0; secondElectron < evt_nels; ++secondElectron) {
      if ( secondElectron == firstElectron ) continue;
      // e
      for (unsigned int thirdElectron = 0; thirdElectron < evt_nels; ++thirdElectron) {
	if ( thirdElectron == firstElectron || thirdElectron == secondElectron) continue;	      
	sorter[0] = firstElectron;
	sorter[1] = secondElectron;
	sorter[2] = thirdElectron;
	std::sort(sorter,sorter+3);
	unsigned int candIndex = encodeTrileptonCandidate(sorter[0],
							  sorter[1],
							  sorter[2]);
	if ( find(cand.begin(),cand.end(), candIndex) == cand.end() ) {
	  cand.push_back(candIndex);
	  // check for charge to identify bucket number
	  int charge = els_charge->at(sorter[0]) + els_charge->at(sorter[1]) + els_charge->at(sorter[2]);
	  if ( charge == -3 ) {
	    // e-e-e-
	    vector_hyp_trilep_bucket->push_back(3);
	  } else if ( charge == -1) {
	    // e+e-e-
	    vector_hyp_trilep_bucket->push_back(2);
	  } else if ( charge == 1 ) {
	    // e+e+e-
	    vector_hyp_trilep_bucket->push_back(1);
	  } else if ( charge == 3 ) {
	    // e+e+e+
	    vector_hyp_trilep_bucket->push_back(0);
	  } else {
	    edm::LogError("HypTrilepMaker") << "combineTriLeptons eee: charge combination could not be identified!!!";
	  }
	  vector_hyp_trilep_first_type->push_back(els_charge->at(sorter[0])*2);
	  vector_hyp_trilep_first_index->push_back(sorter[0]);
	  vector_hyp_trilep_second_type->push_back(els_charge->at(sorter[1])*2);
	  vector_hyp_trilep_second_index->push_back(sorter[1]);
	  vector_hyp_trilep_third_type->push_back(els_charge->at(sorter[2])*2);
	  vector_hyp_trilep_third_index->push_back(sorter[2]);
	}
      }
    }
  }

  // put containers into event
  iEvent.put(vector_hyp_trilep_bucket,"hyptrilepbucket");
  iEvent.put(vector_hyp_trilep_first_type,"hyptrilepfirsttype");
  iEvent.put(vector_hyp_trilep_first_index,"hyptrilepfirstindex");
  iEvent.put(vector_hyp_trilep_second_type,"hyptrilepsecondtype");
  iEvent.put(vector_hyp_trilep_second_index,"hyptrilepsecondindex");
  iEvent.put(vector_hyp_trilep_third_type,"hyptrilepthirdtype");
  iEvent.put(vector_hyp_trilep_third_index,"hyptrilepthirdindex");
}

// ------------ method called once each job just before starting event loop  ------------
void 
HypTrilepMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HypTrilepMaker::endJob() {
}

unsigned int HypTrilepMaker::encodeTrileptonCandidate(unsigned int first, unsigned int second, unsigned int third) {
  // encode trilepton candidate according to
  //
  // trilepton candidate is identified by coded unsigned integer: BBCCDD
  //
  // BB: index of first lepton from els and mus collection
  // CC: index of second lepton from els and mus collection
  // DD: index of third lepton from els and mus collection
  //
  // flavors are ordered: first mu, then el
  // same flavor leptons in candidate are ordered by increasing index
  //

  return first * 10000 + second * 100 + third;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HypTrilepMaker);

