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
m+m+m+: 0
m+m+m-: 1
m+m-m-: 2
m-m-m-: 3
m+m+e+: 4
m+m+e-: 5
m+m-e+: 6
m+m-e-: 7
m-m-e+: 8
m-m-e-: 9
m+e+e+: 10
m+e+e-: 11
m+e-e-: 12
m-e+e+: 13
m-e+e-: 14
m-e-e-: 15
e+e+e+: 16
e+e+e-: 17
e+e-e-: 18
e-e-e-: 19
- trilepton hypothesis first/second/third follows bucket definition. In the case of the same type, ordering follows descending index number within muon and electron collection.

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: HypTrilepMaker.cc,v 1.9 2009/05/21 01:09:14 fgolf Exp $
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
#include "CMS2/NtupleMaker/interface/METUtilities.h"
#include "CMS2/NtupleMaker/interface/JetUtilities.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//
// class decleration
//

//
// constructors and destructor
//
HypTrilepMaker::HypTrilepMaker(const edm::ParameterSet& iConfig)
{
  // parameters from configuration
  muonsInputTag = iConfig.getParameter<edm::InputTag>("muonsInputTag");
  electronsInputTag = iConfig.getParameter<edm::InputTag>("electronsInputTag");
  jetsInputTag = iConfig.getParameter<edm::InputTag>("jetsInputTag");
  trksInputTag = iConfig.getParameter<edm::InputTag>("trksInputTag");
  hypJetMinEtaCut = iConfig.getParameter<double>("hypJetMinEtaCut");
  hypJetMaxEtaCut = iConfig.getParameter<double>("hypJetMaxEtaCut");
  hypJetMinPtCut = iConfig.getParameter<double>("hypJetMinPtCut");
  tightptcut = iConfig.getParameter<double>("TightLepton_PtCut");
  looseptcut = iConfig.getParameter<double>("LooseLepton_PtCut");

  // product of this EDProducer
  // 
  // trilepton hyptothesis
  //
  produces<std::vector<unsigned int> >                ("hyptrilepbucket").setBranchAlias("hyp_trilep_bucket");            // trilepton bucket
  produces<std::vector<int> >                         ("hyptrilepfirsttype").setBranchAlias("hyp_trilep_first_type");     // type of the first lepton in the trilepton hypothesis (1: muon, 2: electron)
  produces<std::vector<unsigned int> >                ("hyptrilepfirstindex").setBranchAlias("hyp_trilep_first_index");   // index of first lepton in lepton collection
  produces<std::vector<int> >                         ("hyptrilepsecondtype").setBranchAlias("hyp_trilep_second_type");   // type of the second lepton in the trilepton hypothesis (1: muon, 2: electron)
  produces<std::vector<unsigned int> >                ("hyptrilepsecondindex").setBranchAlias("hyp_trilep_second_index"); // index of second lepton in lepton collection
  produces<std::vector<int> >                         ("hyptrilepthirdtype").setBranchAlias("hyp_trilep_third_type");     // type of the third lepton in the trilepton hypothesis (1: muon, 2: electron)
  produces<std::vector<unsigned int> >                ("hyptrilepthirdindex").setBranchAlias("hyp_trilep_third_index");   // index of third lepton in lepton collection
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

  // output collections
  std::auto_ptr<std::vector<unsigned int> > vector_hyp_trilep_bucket(new std::vector<unsigned int>);
  std::auto_ptr<std::vector<int> > vector_hyp_trilep_first_type(new std::vector<int>);
  std::auto_ptr<std::vector<unsigned int> > vector_hyp_trilep_first_index(new std::vector<unsigned int>);
  std::auto_ptr<std::vector<int> > vector_hyp_trilep_second_type(new std::vector<int>);
  std::auto_ptr<std::vector<unsigned int> > vector_hyp_trilep_second_index(new std::vector<unsigned int>);
  std::auto_ptr<std::vector<int> > vector_hyp_trilep_third_type(new std::vector<int>);
  std::auto_ptr<std::vector<unsigned int> > vector_hyp_trilep_third_index(new std::vector<unsigned int>);

  // muon charge
  edm::InputTag mus_charge_tag(muonsInputTag.label(),"muscharge");
  edm::Handle<std::vector<int> > mus_charge_h;
  iEvent.getByLabel(mus_charge_tag, mus_charge_h);
  const std::vector<int> *mus_charge = mus_charge_h.product();

  //muon p4
  edm::InputTag mus_p4_tag(muonsInputTag.label(),"musp4");
  edm::Handle<std::vector<LorentzVector> > mus_p4_h;
  iEvent.getByLabel(mus_p4_tag, mus_p4_h);
  const std::vector<LorentzVector> *mus_p4 = mus_p4_h.product();

  // electron charge
  edm::InputTag els_charge_tag(electronsInputTag.label(),"elscharge");
  edm::Handle<std::vector<int> > els_charge_h;
  iEvent.getByLabel(els_charge_tag, els_charge_h);
  const std::vector<int> *els_charge = els_charge_h.product();

  // electron p4
  edm::InputTag els_p4_tag(electronsInputTag.label(),"elsp4");
  edm::Handle<std::vector<LorentzVector> > els_p4_h;
  iEvent.getByLabel(els_p4_tag, els_p4_h);
  const std::vector<LorentzVector> *els_p4 = els_p4_h.product();

  //jet p4
  edm::InputTag jets_p4_tag(jetsInputTag.label(), "jptsp4");
  edm::Handle<std::vector<LorentzVector> > jets_p4_h;
  iEvent.getByLabel(jets_p4_tag, jets_p4_h);
  const std::vector<LorentzVector> *jets_p4 = jets_p4_h.product();


  // number of electrons
  unsigned int evt_nels = els_charge->size();

  // number of electrons
  unsigned int evt_nmus = mus_charge->size();

  // check for numbers of electrons/muons
  // if more than 99 electrons of 99 muons, skip event and fill empty vectors
  if ( evt_nels > 99 ) {
    edm::LogWarning("HypTrilepMaker") << "more than 99 electrons, skipping event!!!";
    // put empty containers into event
    iEvent.put(vector_hyp_trilep_bucket,"hyptrilepbucket");
    iEvent.put(vector_hyp_trilep_first_type,"hyptrilepfirsttype");
    iEvent.put(vector_hyp_trilep_first_index,"hyptrilepfirstindex");
    iEvent.put(vector_hyp_trilep_second_type,"hyptrilepsecondtype");
    iEvent.put(vector_hyp_trilep_second_index,"hyptrilepsecondindex");
    iEvent.put(vector_hyp_trilep_third_type,"hyptrilepthirdtype");
    iEvent.put(vector_hyp_trilep_third_index,"hyptrilepthirdindex");
    return;
  } else if ( evt_nmus > 99 ) {
    edm::LogWarning("HypTrilepMaker") << "more than 99 muons, skipping event!!!";
    // put empty containers into event
    iEvent.put(vector_hyp_trilep_bucket,"hyptrilepbucket");
    iEvent.put(vector_hyp_trilep_first_type,"hyptrilepfirsttype");
    iEvent.put(vector_hyp_trilep_first_index,"hyptrilepfirstindex");
    iEvent.put(vector_hyp_trilep_second_type,"hyptrilepsecondtype");
    iEvent.put(vector_hyp_trilep_second_index,"hyptrilepsecondindex");
    iEvent.put(vector_hyp_trilep_third_type,"hyptrilepthirdtype");
    iEvent.put(vector_hyp_trilep_third_index,"hyptrilepthirdindex");
    return;
  }

  std::vector<unsigned int> cand;

  // array for sorts
  unsigned int sorter[3] = {0,0,0};

  // m
  for (unsigned int firstMuon = 0; firstMuon < evt_nmus; ++firstMuon) {
    // m
    for (unsigned int secondMuon = firstMuon+1; secondMuon < evt_nmus; ++secondMuon) {
      if ( secondMuon == firstMuon ) continue;
      // m
      for (unsigned int thirdMuon = secondMuon+1; thirdMuon < evt_nmus; ++thirdMuon) {
	if ( thirdMuon == firstMuon || thirdMuon == secondMuon ) continue;

	// hyp lepton pt cuts
	// check that all leptons have >= looseptcut
	if ( mus_p4->at(firstMuon).Pt() < looseptcut ||
	     mus_p4->at(secondMuon).Pt() < looseptcut ||
	     mus_p4->at(thirdMuon).Pt() < looseptcut ) continue;
	// check that at least one lepton has >= tightptcut
	if ( mus_p4->at(firstMuon).Pt() < tightptcut &&
	     mus_p4->at(secondMuon).Pt() < tightptcut &&
	     mus_p4->at(thirdMuon).Pt() < tightptcut ) continue;

	sorter[0] = firstMuon;
	sorter[1] = secondMuon;
	sorter[2] = thirdMuon;
	std::sort(sorter,sorter+3);
	unsigned int candIndex = encodeTrileptonCandidate(0,
							  sorter[0],
							  sorter[1],
							  sorter[2]);
	if ( find(cand.begin(),cand.end(), candIndex) == cand.end() ) {
	  cand.push_back(candIndex);
	  int charge = mus_charge->at(sorter[0]) + mus_charge->at(sorter[1]) + mus_charge->at(sorter[2]);
	  if ( charge == 3 ) {
	    // m+m+m+
	    vector_hyp_trilep_bucket->push_back(0);
	  } else if ( charge == 1) {
	    // m+m+m-
	    vector_hyp_trilep_bucket->push_back(1);
	  } else if ( charge == -1 ) {
	    // m+m-m-
	    vector_hyp_trilep_bucket->push_back(2);
	  } else if ( charge == -3 ) {
	    // m-m-m-
	    vector_hyp_trilep_bucket->push_back(3);
	  } else {
	    edm::LogError("HypTrilepMaker") << "combineTriLeptons mmm : charge combination could not be identified!!!";
	  }

	  // store jet indices which pass cuts
	  std::vector<int>  jets_index;
	  for(unsigned int i = 0; i<jets_p4->size(); ++i) {
	    // jet pt cut, if jets corrected, anti-correct cut value as cut is on uncorrected jets
	    double hypJetMinPtCutAdapted = hypJetMinPtCut;

	    if ( jets_p4->at(i).eta() >= hypJetMaxEtaCut ) continue;
	    if ( jets_p4->at(i).eta() <= hypJetMinEtaCut ) continue;
	    if ( jets_p4->at(i).Pt() <= hypJetMinPtCutAdapted ) continue;
	    jets_index.push_back(i);
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

	// hyp lepton pt cuts
	// check that all leptons have >= looseptcut
	if ( mus_p4->at(firstMuon).Pt() < looseptcut ||
	     mus_p4->at(secondMuon).Pt() < looseptcut ||
	     els_p4->at(thirdElectron).Pt() < looseptcut ) continue;
	// check that at least one lepton has >= tightptcut
	if ( mus_p4->at(firstMuon).Pt() < tightptcut &&
	     mus_p4->at(secondMuon).Pt() < tightptcut &&
	     els_p4->at(thirdElectron).Pt() < tightptcut ) continue;

	// order muon indices
	sorter[0] = firstMuon;
	sorter[1] = secondMuon;
	std::sort(sorter,sorter+2);
	// add electron
	sorter[2] = thirdElectron;
	unsigned int candIndex = encodeTrileptonCandidate(1,
							  sorter[0],
							  sorter[1],
							  sorter[2]);
	if ( find(cand.begin(),cand.end(), candIndex) == cand.end() ) {
	  cand.push_back(candIndex);
	  int charge = mus_charge->at(sorter[0]) + mus_charge->at(sorter[1]);
	  if ( charge == 2 ) {
	    if ( els_charge->at(sorter[2]) > 0 ) {
	      // m+m+e+
	      vector_hyp_trilep_bucket->push_back(4);
	    } else {
	      // m+m+e-
	      vector_hyp_trilep_bucket->push_back(5);
	    }
	  } else if ( charge == 0) {
	    if ( els_charge->at(sorter[2]) > 0 ) {
	      // m+m-e+
	      vector_hyp_trilep_bucket->push_back(6);
	    } else {
	      // m+m-e-
	      vector_hyp_trilep_bucket->push_back(7);
	    }
	  } else if ( charge == -2 ) {
	    if ( els_charge->at(sorter[2]) > 0 ) {
	      // m-m-e+
	      vector_hyp_trilep_bucket->push_back(8);
	    } else {
	      // m-m-e-
	      vector_hyp_trilep_bucket->push_back(9);
	    }
	  } else {
	    edm::LogError("HypTrilepMaker") << "combineTriLeptons mme: charge combination could not be identified!!!";
	  }

	  // store jet indices which pass cuts
	  std::vector<int>  jets_index;
	  for(unsigned int i = 0; i<jets_p4->size(); ++i) {
	    // jet pt cut, if jets corrected, anti-correct cut value as cut is on uncorrected jets
	    double hypJetMinPtCutAdapted = hypJetMinPtCut;

	    // jet pre-selection
	    if ( jets_p4->at(i).eta() >= hypJetMaxEtaCut ) continue;
	    if ( jets_p4->at(i).eta() <= hypJetMinEtaCut ) continue;
	    if ( jets_p4->at(i).Pt() <= hypJetMinPtCutAdapted ) continue;

	    // veto electron jets
	    if(!JetUtilities::testJetForElectrons(jets_p4->at(i), els_p4->at(sorter[2]))) continue;

	    jets_index.push_back(i);
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
      for (unsigned int thirdElectron = secondElectron+1; thirdElectron < evt_nels; ++thirdElectron) {
	if ( thirdElectron == secondElectron ) continue;

	// hyp lepton pt cuts
	// check that all leptons have >= looseptcut
	if ( mus_p4->at(firstMuon).Pt() < looseptcut ||
	     els_p4->at(secondElectron).Pt() < looseptcut ||
	     els_p4->at(thirdElectron).Pt() < looseptcut ) continue;
	// check that at least one lepton has >= tightptcut
	if ( mus_p4->at(firstMuon).Pt() < tightptcut &&
	     els_p4->at(secondElectron).Pt() < tightptcut &&
	     els_p4->at(thirdElectron).Pt() < tightptcut ) continue;


	// order electron indices
	sorter[0] = secondElectron;
	sorter[1] = thirdElectron;
	std::sort(sorter,sorter+2);
	// add muon in first place
	sorter[2] = sorter[1];
	sorter[1] = sorter[0];
	sorter[0] = firstMuon;
	unsigned int candIndex = encodeTrileptonCandidate(2,
							  sorter[0],
							  sorter[1],
							  sorter[2]);
	if ( find(cand.begin(),cand.end(), candIndex) == cand.end() ) {
	  cand.push_back(candIndex);
	  // check for charge to identify bucket number
	  int charge = els_charge->at(sorter[1]) + els_charge->at(sorter[2]);
	  if ( mus_charge->at(sorter[0]) > 0 ) {
	    if ( charge == 2 ) {
	      // m+e+e+
	      vector_hyp_trilep_bucket->push_back(10);
	    } else if ( charge == 0) {
	      // m+e+e-
	      vector_hyp_trilep_bucket->push_back(11);
	    } else if ( charge == -2 ) {
	      // m+e-e-
	      vector_hyp_trilep_bucket->push_back(12);
	    } else {
	      edm::LogError("HypTrilepMaker") << "combineTriLeptons m+ee: charge combination could not be identified!!!";
	    }
	  } else if ( mus_charge->at(sorter[0]) < 0 ) {
	    if ( charge == 2 ) {
	      // m-e+e+
	      vector_hyp_trilep_bucket->push_back(13);
	    } else if ( charge == 0) {
	      // m-e+e-
	      vector_hyp_trilep_bucket->push_back(14);
	    } else if ( charge == -2 ) {
	      // m-e-e-
	      vector_hyp_trilep_bucket->push_back(15);
	    } else {
	      edm::LogError("HypTrilepMaker") << "combineTriLeptons m-ee: charge combination could not be identified!!!";
	    }
	  } else {
	    edm::LogError("HypTrilepMaker") << "combineTriLeptons mee: charge combination could not be identified!!!";
	  }
	  }

	  // store jet indices which pass cuts
	  std::vector<int>  jets_index;
	  for(unsigned int i = 0; i<jets_p4->size(); ++i) {
	    // jet pt cut, if jets corrected, anti-correct cut value as cut is on uncorrected jets
	    double hypJetMinPtCutAdapted = hypJetMinPtCut;

	    // jet pre-selection
	    if ( jets_p4->at(i).eta() >= hypJetMaxEtaCut ) continue;
	    if ( jets_p4->at(i).eta() <= hypJetMinEtaCut ) continue;
	    if ( jets_p4->at(i).Pt() <= hypJetMinPtCutAdapted ) continue;

	    // veto electron jets
	    if(!JetUtilities::testJetForElectrons(jets_p4->at(i), els_p4->at(sorter[1]))) continue;
	    if(!JetUtilities::testJetForElectrons(jets_p4->at(i), els_p4->at(sorter[2]))) continue;

	    jets_index.push_back(i);
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

  
  // e
for (unsigned int firstElectron = 0; firstElectron < evt_nels; ++firstElectron) {
    // e
    for (unsigned int secondElectron = firstElectron+1; secondElectron < evt_nels; ++secondElectron) {
      if ( secondElectron == firstElectron ) continue;
      // e
      for (unsigned int thirdElectron = secondElectron+1; thirdElectron < evt_nels; ++thirdElectron) {
	if ( thirdElectron == firstElectron || thirdElectron == secondElectron) continue;	      

	// hyp lepton pt cuts
	// check that all leptons have >= looseptcut
	if ( els_p4->at(firstElectron).Pt() < looseptcut ||
	     els_p4->at(secondElectron).Pt() < looseptcut ||
	     els_p4->at(thirdElectron).Pt() < looseptcut ) continue;
	// check that at least one lepton has >= tightptcut
	if ( els_p4->at(firstElectron).Pt() < tightptcut &&
	     els_p4->at(secondElectron).Pt() < tightptcut &&
	     els_p4->at(thirdElectron).Pt() < tightptcut ) continue;

	sorter[0] = firstElectron;
	sorter[1] = secondElectron;
	sorter[2] = thirdElectron;
	std::sort(sorter,sorter+3);
	unsigned int candIndex = encodeTrileptonCandidate(3,
							  sorter[0],
							  sorter[1],
							  sorter[2]);
	if ( find(cand.begin(),cand.end(), candIndex) == cand.end() ) {
	  cand.push_back(candIndex);
	  // check for charge to identify bucket number
	  int charge = els_charge->at(sorter[0]) + els_charge->at(sorter[1]) + els_charge->at(sorter[2]);
	  if ( charge == 3 ) {
	    // e+e+e+
	    vector_hyp_trilep_bucket->push_back(16);
	  } else if ( charge == 1) {
	    // e+e+e-
	    vector_hyp_trilep_bucket->push_back(17);
	  } else if ( charge == -1 ) {
	    // e+e-e-
	    vector_hyp_trilep_bucket->push_back(18);
	  } else if ( charge == -3 ) {
	    // e-e-e-
	    vector_hyp_trilep_bucket->push_back(19);
	  } else {
	    edm::LogError("HypTrilepMaker") << "combineTriLeptons eee: charge combination could not be identified!!!";
	  }

	  // store jet indices which pass cuts
	  std::vector<int>  jets_index;
	  for(unsigned int i = 0; i<jets_p4->size(); ++i) {
	    // jet pt cut, if jets corrected, anti-correct cut value as cut is on uncorrected jets
	    double hypJetMinPtCutAdapted = hypJetMinPtCut;

	    // jet pre-selection
	    if ( jets_p4->at(i).eta() >= hypJetMaxEtaCut ) continue;
	    if ( jets_p4->at(i).eta() <= hypJetMinEtaCut ) continue;
	    if ( jets_p4->at(i).Pt() <= hypJetMinPtCutAdapted ) continue;

	    // veto electron jets
	    if(!JetUtilities::testJetForElectrons(jets_p4->at(i), els_p4->at(sorter[0]))) continue;
	    if(!JetUtilities::testJetForElectrons(jets_p4->at(i), els_p4->at(sorter[1]))) continue;
	    if(!JetUtilities::testJetForElectrons(jets_p4->at(i), els_p4->at(sorter[2]))) continue;

	    jets_index.push_back(i);
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

unsigned int HypTrilepMaker::encodeTrileptonCandidate(unsigned int combination, unsigned int first, unsigned int second, unsigned int third) {
  // encode trilepton candidate according to
  //
  // trilepton candidate is identified by coded unsigned integer: ABBCCDD
  //
  //  A: trilepton combination: mmm:0, mme: 1, mee: 2, eee: 3
  // BB: index of first lepton from els and mus collection
  // CC: index of second lepton from els and mus collection
  // DD: index of third lepton from els and mus collection
  //
  // flavors are ordered: first mu, then el
  // same flavor leptons in candidate are ordered by increasing index
  //

  return combination * 1000000 + first * 10000 + second * 100 + third;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HypTrilepMaker);

