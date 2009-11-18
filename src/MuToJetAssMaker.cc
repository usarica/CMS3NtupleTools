// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuToJetAssMaker
// 
/**\class MuToJetAssMaker MuToJetAssMaker.cc CMS2/NtupleMaker/src/MuToJetAssMaker.cc

 Description: make associations between jets and muons

*/
//
// Original Author:  Frank Golf
//         Created:  Wed Jun 25 18:32:24 UTC 2008
// $Id: MuToJetAssMaker.cc,v 1.4 2009/11/18 21:46:25 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/MuToJetAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

MuToJetAssMaker::MuToJetAssMaker(const edm::ParameterSet& iConfig)
     : m_minDR(iConfig.getParameter<double>("minDR"))
{
     produces<vector<int>   >("musclosestJet").setBranchAlias("mus_closestJet");	// muon closest to jet
     produces<vector<float> >("musjetdr"     ).setBranchAlias("mus_jetdr"     );     
}

void MuToJetAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>   > vector_mus_closestJet(new vector<int>  );
     std::auto_ptr<vector<float> > vector_mus_jetdr     (new vector<float>);

     // get muons
     Handle<vector<LorentzVector> > mus_p4_h;
     iEvent.getByLabel("muonMaker", "musp4", mus_p4_h);  

     // get jet p4's
     Handle<vector<LorentzVector> > jets_p4_h;
     iEvent.getByLabel("jetMaker", "jetsp4", jets_p4_h);
     
     // loop over all muons
     for(vector<LorentzVector>::const_iterator mus_it = mus_p4_h->begin(),
	 mus_end = mus_p4_h->end();
	 mus_it != mus_end; ++mus_it) {
	 
       double mu_eta = mus_it->Eta();
       double mu_phi = mus_it->Phi();
       
       double minDR = m_minDR;
       unsigned int i = 0;
       int index      = -9999; 

       for(vector<LorentzVector>::const_iterator jets_it = jets_p4_h->begin(),
	   jets_end = jets_p4_h->end();
	   jets_it != jets_end; ++jets_it, ++i) {

	 double jet_eta = jets_it->Eta();
	 double jet_phi = jets_it->Phi();

	 double dR = deltaR(mu_eta, mu_phi, jet_eta, jet_phi);
	 if(dR < minDR) {
	   minDR = dR;
	   index = i;
	 }
       }

       // fill vector
       vector_mus_closestJet->push_back(index);
       vector_mus_jetdr     ->push_back(minDR);
     }

     // store vectors
     iEvent.put(vector_mus_closestJet, "musclosestJet");
     iEvent.put(vector_mus_jetdr     , "musjetdr"     );
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuToJetAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuToJetAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuToJetAssMaker);
