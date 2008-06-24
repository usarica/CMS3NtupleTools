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
// $Id: HypTrilepMaker.cc,v 1.1 2008/06/24 00:39:34 gutsche Exp $
//
//

// system include files
#include <memory>
#include <vector>
#include <algorithm>

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
  jetsInputTag = iConfig.getParameter<edm::InputTag>("jetsInputTag");
  muonsInputTag = iConfig.getParameter<edm::InputTag>("muonsInputTag");
  electronsInputTag = iConfig.getParameter<edm::InputTag>("electronsInputTag");
  elToMuAssInputTag = iConfig.getParameter<edm::InputTag>("elToMuAssInputTag");
  metInputTag = iConfig.getParameter<edm::InputTag>("metInputTag");

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

  // muon p4
  edm::InputTag mus_p4_tag(muonsInputTag.label(),"musp4");
  edm::Handle<std::vector<LorentzVector> > mus_p4_h;
  iEvent.getByLabel(mus_p4_tag, mus_p4_h);
  const std::vector<LorentzVector> *mus_p4 = mus_p4_h.product();

  // muon charge
  edm::InputTag mus_charge_tag(muonsInputTag.label(),"muscharge");
  edm::Handle<std::vector<int> > mus_charge_h;
  iEvent.getByLabel(mus_charge_tag, mus_charge_h);
  const std::vector<int> *mus_charge = mus_charge_h.product();

  // muon gfit chi2
  edm::InputTag mus_gfit_chi2_tag(muonsInputTag.label(),"musgfitchi2");
  edm::Handle<std::vector<float> > mus_gfit_chi2_h;
  iEvent.getByLabel(mus_gfit_chi2_tag, mus_gfit_chi2_h);
  const std::vector<float> *mus_gfit_chi2 = mus_gfit_chi2_h.product();

  // muon gfit ndof
  edm::InputTag mus_gfit_ndof_tag(muonsInputTag.label(),"musgfitndof");
  edm::Handle<std::vector<float> > mus_gfit_ndof_h;
  iEvent.getByLabel(mus_gfit_ndof_tag, mus_gfit_ndof_h);
  const std::vector<float> *mus_gfit_ndof = mus_gfit_ndof_h.product();

  // muon d0
  edm::InputTag mus_d0_tag(muonsInputTag.label(),"musd0");
  edm::Handle<std::vector<float> > mus_d0_h;
  iEvent.getByLabel(mus_d0_tag, mus_d0_h);
  const std::vector<float> *mus_d0 = mus_d0_h.product();

  // muon validHits
  edm::InputTag mus_validHits_tag(muonsInputTag.label(),"musvalidHits");
  edm::Handle<std::vector<int> > mus_validHits_h;
  iEvent.getByLabel(mus_validHits_tag, mus_validHits_h);
  const std::vector<int> *mus_validHits = mus_validHits_h.product();

  // muon iso03 sumPt
  edm::InputTag mus_iso03_sumPt_tag(muonsInputTag.label(),"musiso03sumPt");
  edm::Handle<std::vector<float> > mus_iso03_sumPt_h;
  iEvent.getByLabel(mus_iso03_sumPt_tag, mus_iso03_sumPt_h);
  const std::vector<float> *mus_iso03_sumPt = mus_iso03_sumPt_h.product();

  // muon iso03 emEt
  edm::InputTag mus_iso03_emEt_tag(muonsInputTag.label(),"musiso03emEt");
  edm::Handle<std::vector<float> > mus_iso03_emEt_h;
  iEvent.getByLabel(mus_iso03_emEt_tag, mus_iso03_emEt_h);
  const std::vector<float> *mus_iso03_emEt = mus_iso03_emEt_h.product();

  // muon iso03 hadEt
  edm::InputTag mus_iso03_hadEt_tag(muonsInputTag.label(),"musiso03hadEt");
  edm::Handle<std::vector<float> > mus_iso03_hadEt_h;
  iEvent.getByLabel(mus_iso03_hadEt_tag, mus_iso03_hadEt_h);
  const std::vector<float> *mus_iso03_hadEt = mus_iso03_hadEt_h.product();

  // electron p4
  edm::InputTag els_p4_tag(electronsInputTag.label(),"elsp4");
  edm::Handle<std::vector<LorentzVector> > els_p4_h;
  iEvent.getByLabel(els_p4_tag, els_p4_h);
  const std::vector<LorentzVector> *els_p4 = els_p4_h.product();

  // electron charge
  edm::InputTag els_charge_tag(electronsInputTag.label(),"elscharge");
  edm::Handle<std::vector<int> > els_charge_h;
  iEvent.getByLabel(els_charge_tag, els_charge_h);
  const std::vector<int> *els_charge = els_charge_h.product();

  // electron tightId
  edm::InputTag els_tightId_tag(electronsInputTag.label(),"elstightId");
  edm::Handle<std::vector<int> > els_tightId_h;
  iEvent.getByLabel(els_tightId_tag, els_tightId_h);
  const std::vector<int> *els_tightId = els_tightId_h.product();

  // electron d0
  edm::InputTag els_d0_tag(electronsInputTag.label(),"elsd0");
  edm::Handle<std::vector<float> > els_d0_h;
  iEvent.getByLabel(els_d0_tag, els_d0_h);
  const std::vector<float> *els_d0 = els_d0_h.product();

  // electron tkIso
  edm::InputTag els_tkIso_tag(electronsInputTag.label(),"elstkIso");
  edm::Handle<std::vector<float> > els_tkIso_h;
  iEvent.getByLabel(els_tkIso_tag, els_tkIso_h);
  const std::vector<float> *els_tkIso = els_tkIso_h.product();

  // electron closestMuon
  edm::InputTag els_closestMuon_tag(elToMuAssInputTag.label(),"elsclosestMuon");
  edm::Handle<std::vector<int> > els_closestMuon_h;
  iEvent.getByLabel(els_closestMuon_tag, els_closestMuon_h);
  const std::vector<int> *els_closestMuon = els_closestMuon_h.product();

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
  unsigned int evt_nels = els_p4->size();

  // number of electrons
  unsigned int evt_nmus = mus_p4->size();

  // m
  for (unsigned int firstMuon = 0; firstMuon < evt_nmus; ++firstMuon) {
    if ( goodMuonIsolated(mus_gfit_chi2->at(firstMuon), 
			  mus_gfit_ndof->at(firstMuon), 
			  mus_d0->at(firstMuon), 
			  mus_validHits->at(firstMuon), 
			  mus_iso03_sumPt->at(firstMuon), 
			  mus_iso03_emEt->at(firstMuon), 
			  mus_iso03_hadEt->at(firstMuon), 
			  mus_p4->at(firstMuon)) ) {
      // m
      for (unsigned int secondMuon = 0; secondMuon < evt_nmus; ++secondMuon) {
	if ( secondMuon == firstMuon ) continue;
	if ( goodMuonIsolated(mus_gfit_chi2->at(secondMuon), 
			      mus_gfit_ndof->at(secondMuon), 
			      mus_d0->at(secondMuon), 
			      mus_validHits->at(secondMuon), 
			      mus_iso03_sumPt->at(secondMuon), 
			      mus_iso03_emEt->at(secondMuon), 
			      mus_iso03_hadEt->at(secondMuon), 
			      mus_p4->at(secondMuon)) ) {
	  // m
	  for (unsigned int thirdMuon = 0; thirdMuon < evt_nmus; ++thirdMuon) {
	    if ( thirdMuon == firstMuon || thirdMuon == secondMuon ) continue;
	    if ( goodMuonIsolated(mus_gfit_chi2->at(thirdMuon), 
				  mus_gfit_ndof->at(thirdMuon), 
				  mus_d0->at(thirdMuon), 
				  mus_validHits->at(thirdMuon), 
				  mus_iso03_sumPt->at(thirdMuon), 
				  mus_iso03_emEt->at(thirdMuon), 
				  mus_iso03_hadEt->at(thirdMuon), 
				  mus_p4->at(thirdMuon)) ) {
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
	  }
	  // e
	  for (unsigned int thirdElectron = 0; thirdElectron < evt_nels; ++thirdElectron) {
	    if ( goodElectronIsolated(els_tightId->at(thirdElectron), 
				      els_closestMuon->at(thirdElectron), 
				      els_d0->at(thirdElectron), 
				      els_tkIso->at(thirdElectron), 
				      els_p4->at(thirdElectron)) ) {
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
	}
      }
      // e
      for (unsigned int secondElectron = 0; secondElectron < evt_nels; ++secondElectron) {
	if ( goodElectronIsolated(els_tightId->at(secondElectron), 
				  els_closestMuon->at(secondElectron), 
				  els_d0->at(secondElectron), 
				  els_tkIso->at(secondElectron), 
				  els_p4->at(secondElectron)) ) {
	  // e
	  for (unsigned int thirdElectron = 0; thirdElectron < evt_nels; ++thirdElectron) {
	    if ( thirdElectron == secondElectron ) continue;
	    if ( goodElectronIsolated(els_tightId->at(thirdElectron), 
				      els_closestMuon->at(thirdElectron), 
				      els_d0->at(thirdElectron), 
				      els_tkIso->at(thirdElectron), 
				      els_p4->at(thirdElectron)) ) {
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
		vector_hyp_trilep_second_type->push_back(mus_charge->at(sorter[1])*2);
		vector_hyp_trilep_second_index->push_back(sorter[1]);
		vector_hyp_trilep_third_type->push_back(els_charge->at(sorter[2])*2);
		vector_hyp_trilep_third_index->push_back(sorter[2]);
               

	      }
	    }
	  }
	}
      }
    }
  }

  // e
  for (unsigned int firstElectron = 0; firstElectron < evt_nels; ++firstElectron) {
    if ( goodElectronIsolated(els_tightId->at(firstElectron), 
			      els_closestMuon->at(firstElectron), 
			      els_d0->at(firstElectron), 
			      els_tkIso->at(firstElectron), 
			      els_p4->at(firstElectron)) ) {
      // e
      for (unsigned int secondElectron = 0; secondElectron < evt_nels; ++secondElectron) {
	if ( secondElectron == firstElectron ) continue;
	if ( goodElectronIsolated(els_tightId->at(secondElectron), 
				  els_closestMuon->at(secondElectron), 
				  els_d0->at(secondElectron), 
				  els_tkIso->at(secondElectron), 
				  els_p4->at(secondElectron)) ) {
	  // e
	  for (unsigned int thirdElectron = 0; thirdElectron < evt_nels; ++thirdElectron) {
	    if ( thirdElectron == firstElectron || thirdElectron == secondElectron) continue;
	    if ( goodElectronIsolated(els_tightId->at(thirdElectron), 
				      els_closestMuon->at(thirdElectron), 
				      els_d0->at(thirdElectron), 
				      els_tkIso->at(thirdElectron), 
				      els_p4->at(thirdElectron)) ) {
	      
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
		vector_hyp_trilep_first_type->push_back(mus_charge->at(sorter[0])*2);
		vector_hyp_trilep_first_index->push_back(sorter[0]);
		vector_hyp_trilep_second_type->push_back(mus_charge->at(sorter[1])*2);
		vector_hyp_trilep_second_index->push_back(sorter[1]);
		vector_hyp_trilep_third_type->push_back(els_charge->at(sorter[2])*2);
		vector_hyp_trilep_third_index->push_back(sorter[2]);
	      }
	    }
	  }
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

bool HypTrilepMaker::goodElectronWithoutIsolation(const int els_tightId, 
						  const int els_closestMuon, 
						  const float els_d0) {
  //----------------------------------------------------------------
  // Electron ID without isolation
  //----------------------------------------------------------------
  if ( els_tightId     !=  1)     return false;
  if ( els_closestMuon != -1)     return false;
  if ( std::abs(els_d0)     > 0.025)   return false;
  return true;
}

bool HypTrilepMaker::goodMuonWithoutIsolation(const float mus_gfit_chi2, 
					      const float mus_gfit_ndof, 
					      const float mus_d0, 
					      const int mus_validHits) {
  //----------------------------------------------------------------
  // Muon ID without isolation
  //---------------------------------------------------------------
  if (mus_gfit_chi2/mus_gfit_ndof > 5.)   return false;
  if (std::abs(mus_d0)                 > 0.25) return false;
  if (mus_validHits               < 7)    return false;
  return true;
}

bool HypTrilepMaker::passElectronIsolation(const float els_tkIso, 
					   const LorentzVector &els_p4) {
  //-----------------------------------------------------------
  // Electron Isolation
  //-----------------------------------------------------------
  double sum = els_tkIso;
  double pt  = els_p4.pt();
  if ( pt/(pt+sum) < 0.92) return false;
  return true;  
} 

bool HypTrilepMaker::passMuonIsolation(const float mus_iso03_sumPt, 
				       const float mus_iso03_emEt, 
				       const float mus_iso03_hadEt, 
				       const LorentzVector &mus_p4) {
  //-----------------------------------------------------------
  // Muon Isolation
  //-----------------------------------------------------------
  double sum =  mus_iso03_sumPt + mus_iso03_emEt + mus_iso03_hadEt;
  double pt  = mus_p4.pt(); 
  if ( pt/(pt+sum) < 0.92) return false;
  return true;  
}

bool HypTrilepMaker::goodMuonIsolated(const float mus_gfit_chi2, 
				      const float mus_gfit_ndof, 
				      const float mus_d0, 
				      const int mus_validHits, 
				      const float mus_iso03_sumPt, 
				      const float mus_iso03_emEt, 
				      const float mus_iso03_hadEt, 
				      const LorentzVector &mus_p4) {
  //--------------------------------------------
  // Muon ID with isolation
  //--------------------------------------------
  if (!goodMuonWithoutIsolation(mus_gfit_chi2, 
				mus_gfit_ndof, 
				mus_d0, 
				mus_validHits)) return false;
  if (!passMuonIsolation(mus_iso03_sumPt, 
			 mus_iso03_emEt, 
			 mus_iso03_hadEt, 
			 mus_p4))   return false;
  return true;
}

bool HypTrilepMaker::goodElectronIsolated(const int els_tightId, 
					  const int els_closestMuon, 
					  const float els_d0, 
					  const float els_tkIso, 
					  const LorentzVector &els_p4) {
  //--------------------------------------------
  // Electron ID with isolation
  //--------------------------------------------
  if (!goodElectronWithoutIsolation(els_tightId, 
				    els_closestMuon, 
				    els_d0)) return false;
  if (!passElectronIsolation(els_tkIso, 
			     els_p4))                           return false;
  return true;
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

