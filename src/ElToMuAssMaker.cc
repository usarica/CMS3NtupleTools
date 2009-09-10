// -*- C++ -*-
//
// Package:    ElToMuAssMaker
// Class:      ElToMuAssMaker
// 
/**\class ElToMuAssMaker ElToMuAssMaker.cc CMS2/NtupleMaker/src/ElToMuAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElToMuAssMaker.cc,v 1.7 2009/09/10 10:51:43 fgolf Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/ElToMuAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

ElToMuAssMaker::ElToMuAssMaker(const edm::ParameterSet& iConfig)
     : m_minDR(iConfig.getParameter<double>("minDR"))
{
     produces<vector<int>   >("elsclosestMuon").setBranchAlias("els_closestMuon");	// track index matched to muon
     produces<vector<float> >("elsmusdr"      ).setBranchAlias("els_musdr"      );
}

void ElToMuAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>    > vector_els_closestMuon(new vector<int>   );
     std::auto_ptr<vector<float> > vector_els_musdr      (new vector<float>);
     // get muons
     Handle<vector<LorentzVector> > mus_p4_h;
     iEvent.getByLabel("muonMaker", "musp4", mus_p4_h);  
     //get the muon type
     Handle<vector<int> > mus_type_h;
     iEvent.getByLabel("muonMaker", "mustype", mus_type_h);
     // get track p4's
     Handle<vector<LorentzVector> > els_p4_h;
     iEvent.getByLabel("electronMaker", "elsp4", els_p4_h);

     
     //loop over electrons and find the closest 
     //muon
     for(vector<LorentzVector>::const_iterator els_it = els_p4_h->begin();
	 els_it != els_p4_h->end(); els_it++) {
       
       double el_eta = els_it->Eta();
       double el_phi = els_it->Phi();
       
       double minDR = m_minDR;
       unsigned int i = 0;
       int index = -1; 

       for(vector<LorentzVector>::const_iterator mus_it = mus_p4_h->begin();
	   mus_it != mus_p4_h->end(); mus_it++, i++) {
	 
	 //don't consider this muon if its not a global muon and 
	 //or a tracker muon
	 if(mus_type_h->at(i) == 8)
	   continue;

	 double mu_eta = mus_it->Eta();
	 double mu_phi = mus_it->Phi();

	 double dR = deltaR(el_eta, el_phi, mu_eta, mu_phi);

	 if(dR < minDR) {
	   minDR = dR;
	   index = i;
	 }
       }
       // fill vector
       vector_els_closestMuon->push_back(index);
       vector_els_musdr      ->push_back(minDR);
     }
     // store vectors
     iEvent.put(vector_els_closestMuon, "elsclosestMuon");
     iEvent.put(vector_els_musdr      , "elsmusdr"      );
}

// ------------ method called once each job just before starting event loop  ------------
void 
ElToMuAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElToMuAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElToMuAssMaker);
