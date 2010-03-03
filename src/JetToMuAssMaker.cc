// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      JetToMuAssMaker
// 
/**\class JetToMuAssMaker JetToMuAssMaker.cc CMS2/NtupleMaker/src/JetToMuAssMaker.cc

 Description: make associations between jets and muons

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Tue Jun 17 20:40:42 UTC 2008
// $Id: JetToMuAssMaker.cc,v 1.6 2010/03/03 04:23:54 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/JetToMuAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

JetToMuAssMaker::JetToMuAssMaker(const edm::ParameterSet& iConfig) {
     produces<vector<int> >("jetsclosestMuon").setBranchAlias("jets_closestMuon");	// muon closest to jet
     produces<vector<double> >("jetsclosestMuonDR").setBranchAlias("jets_closestMuon_DR");	// Delta R of muon closest to jet

     m_minDR_       = iConfig.getParameter<double>("minDR"                 );
     jetsInputTag_  = iConfig.getParameter<edm::InputTag>("jetsInputTag"   );
     musInputTag_   = iConfig.getParameter<edm::InputTag>("musInputTag"    );
     
}

void JetToMuAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int> > vector_jets_closestMuon(new vector<int>);        
     std::auto_ptr<vector<double> > vector_jets_closestMuon_DR(new vector<double>);        
     // get muons
     Handle<vector<LorentzVector> > mus_p4_h;
     iEvent.getByLabel(musInputTag_, mus_p4_h);  
     // get track p4's
     Handle<vector<LorentzVector> > jets_p4_h;
     iEvent.getByLabel(jetsInputTag_, jets_p4_h);

     
     // loop over all jets
     vector<LorentzVector>::const_iterator jets_it_end = jets_p4_h->end();
     for(vector<LorentzVector>::const_iterator jets_it = jets_p4_h->begin();
	 jets_it != jets_it_end; jets_it++) {
       
       double jet_eta = jets_it->Eta();
       double jet_phi = jets_it->Phi();
       
       double minDR = m_minDR_;
       unsigned int i = 0;
       int index = -1; 
       for(vector<LorentzVector>::const_iterator mus_it = mus_p4_h->begin();
	   mus_it != mus_p4_h->end(); mus_it++, i++) {
	 
	 double mu_eta = mus_it->Eta();
	 double mu_phi = mus_it->Phi();
	 double dR = deltaR(jet_eta, jet_phi, mu_eta, mu_phi);
	 if(dR < minDR) {
	   minDR = dR;
	   index = i;
	 }
       }
       // fill vector
       vector_jets_closestMuon->push_back(index);
       vector_jets_closestMuon_DR->push_back(minDR);
     }
     // store vectors
     iEvent.put(vector_jets_closestMuon, "jetsclosestMuon");
     iEvent.put(vector_jets_closestMuon_DR, "jetsclosestMuonDR");
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetToMuAssMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetToMuAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetToMuAssMaker);
