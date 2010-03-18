// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ElToJetAssMaker
// 
/**\class ElToJetAssMaker ElToJetAssMaker.cc CMS2/NtupleMaker/src/ElToJetAssMaker.cc

 Description: make associations between jets and electrons

*/
//
// Original Author:  Frank Golf
//         Created:  Wed Jun 25 18:32:24 UTC 2008
// $Id: ElToJetAssMaker.cc,v 1.6 2010/03/18 02:12:02 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/ElToJetAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

ElToJetAssMaker::ElToJetAssMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  
     produces<vector<int>   >(branchprefix+"closestJet").setBranchAlias(aliasprefix_+"_closestJet");	// electron closest to jet
     produces<vector<float> >(branchprefix+"jetdr"     ).setBranchAlias(aliasprefix_+"_jetdr"     );     
     
     elsInputTag_  = iConfig.getParameter<edm::InputTag>("elsInputTag");
     jetsInputTag_ = iConfig.getParameter<edm::InputTag>("jetsInputTag");
     m_minDR_       = iConfig.getParameter<double>       ("minDR");
}

void ElToJetAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>   > vector_els_closestJet(new vector<int>  );
     std::auto_ptr<vector<float> > vector_els_jetdr     (new vector<float>);

     // get electrons
     Handle<vector<LorentzVector> > els_p4_h;
     iEvent.getByLabel(elsInputTag_.label(), "elsp4", els_p4_h);  

     // get jet p4's
     Handle<vector<LorentzVector> > jets_p4_h;
     iEvent.getByLabel(jetsInputTag_.label(), "jetsp4", jets_p4_h);
     
     // loop over all electrons
     for(vector<LorentzVector>::const_iterator els_it = els_p4_h->begin(),
	 els_end = els_p4_h->end();
	 els_it != els_end; ++els_it) {
	 
       double el_eta = els_it->Eta();
       double el_phi = els_it->Phi();
       
       double minDR = m_minDR_;
       unsigned int i = 0;
       int index      = -9999; 

       for(vector<LorentzVector>::const_iterator jets_it = jets_p4_h->begin(),
	   jets_end = jets_p4_h->end();
	   jets_it != jets_end; ++jets_it, ++i) {

	 double jet_eta = jets_it->Eta();
	 double jet_phi = jets_it->Phi();

	 double dR = deltaR(el_eta, el_phi, jet_eta, jet_phi);
	 if(dR < minDR) {
	   minDR = dR;
	   index = i;
	 }
       }

       // fill vector
       vector_els_closestJet->push_back(index);
       vector_els_jetdr     ->push_back(minDR);
     }

     // store vectors
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     iEvent.put(vector_els_closestJet, branchprefix+"closestJet");
     iEvent.put(vector_els_jetdr     , branchprefix+"jetdr"     );
}

// ------------ method called once each job just before starting event loop  ------------
void ElToJetAssMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void ElToJetAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElToJetAssMaker);
