// -*- C++ -*-
//
// Package:    JetToElAssMaker
// Class:      JetToElAssMaker
// 
/**\class JetToElAssMaker JetToElAssMaker.cc CMS2/NtupleMaker/src/JetToElAssMaker.cc

 Description: make associations between jets and electrons

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Tue Jun 17 20:40:42 UTC 2008
// $Id: JetToElAssMaker.cc,v 1.2 2008/09/13 08:07:23 jmuelmen Exp $
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
#include "CMS2/NtupleMaker/interface/JetToElAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

JetToElAssMaker::JetToElAssMaker(const edm::ParameterSet& iConfig)
     : m_minDR(iConfig.getParameter<double>("minDR"))
{
     produces<vector<int> >("jetsclosestElectron").setBranchAlias("jets_closestElectron");		// electron closest to jet
     produces<vector<double> >("jetsclosestElectronDR").setBranchAlias("jets_closestElectron_DR");	// Delta R of electron closest to jet
}

void JetToElAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int> > vector_jets_closestElectron(new vector<int>);        
     std::auto_ptr<vector<double> > vector_jets_closestElectron_DR(new vector<double>);        
     // get muons
     Handle<vector<LorentzVector> > els_p4_h;
     iEvent.getByLabel("electronMaker", "elsp4", els_p4_h);  
     // get track p4's
     Handle<vector<LorentzVector> > jets_p4_h;
     iEvent.getByLabel("jetMaker", "jetsp4", jets_p4_h);

     
     // loop over all jets
     vector<LorentzVector>::const_iterator jets_it_end = jets_p4_h->end();
     for(vector<LorentzVector>::const_iterator jets_it = jets_p4_h->begin();
	 jets_it != jets_it_end; jets_it++) {
       
       double jet_eta = jets_it->Eta();
       double jet_phi = jets_it->Phi();
       
       double minDR = m_minDR;
       unsigned int i = 0;
       int index = -1; 
       for(vector<LorentzVector>::const_iterator els_it = els_p4_h->begin();
	   els_it != els_p4_h->end(); els_it++, i++) {
	 
	 double mu_eta = els_it->Eta();
	 double mu_phi = els_it->Phi();
	 double dR = deltaR(jet_eta, jet_phi, mu_eta, mu_phi);
	 if(dR < minDR) {
	   minDR = dR;
	   index = i;
	 }
       }
       // fill vector
       vector_jets_closestElectron->push_back(index);
       vector_jets_closestElectron_DR->push_back(minDR);
     }
     // store vectors
     iEvent.put(vector_jets_closestElectron, "jetsclosestElectron");
     iEvent.put(vector_jets_closestElectron_DR, "jetsclosestElectronDR");
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetToElAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetToElAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetToElAssMaker);
