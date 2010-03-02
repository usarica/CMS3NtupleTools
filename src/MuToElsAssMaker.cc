// -*- C++ -*-
//
// Package:    MuToElsAssMaker
// Class:      MuToElsAssMaker
// 
/**\class MuToElsAssMaker MuToElsAssMaker.cc CMS2/NtupleMaker/src/MuToElsAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuToElsAssMaker.cc,v 1.5 2010/03/02 19:36:08 fgolf Exp $
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
#include "CMS2/NtupleMaker/interface/MuToElsAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

MuToElsAssMaker::MuToElsAssMaker(const edm::ParameterSet& iConfig)
     : m_minDR(iConfig.getParameter<double>("minDR"))
{
     produces<vector<int>   >("musclosestEle").setBranchAlias("mus_closestEle");	// muon matched to electron
     produces<vector<float> >("museledr"     ).setBranchAlias("mus_eledr"     );
}

void MuToElsAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>   > vector_mus_closestEle(new vector<int>  );
     std::auto_ptr<vector<float> > vector_mus_eledr     (new vector<float>);

     // get muons
     Handle<vector<LorentzVector> > mus_p4_h;
     iEvent.getByLabel("muonMaker", "musp4", mus_p4_h);  

     // get electrons
     Handle<vector<LorentzVector> > els_p4_h;
     iEvent.getByLabel("electronMaker", "elsp4", els_p4_h);
     
     //loop over muons and find the closest electron
     for(vector<LorentzVector>::const_iterator mus_it = mus_p4_h->begin(),
	 mus_end = mus_p4_h->end();
	 mus_it != mus_end; ++mus_it) {
       
       double mu_eta = mus_it->Eta();
       double mu_phi = mus_it->Phi();
       
       double minDR = m_minDR;
       unsigned int i = 0;
       int index      = -9999;
 
       for(vector<LorentzVector>::const_iterator els_it = els_p4_h->begin(),
	   els_end = els_p4_h->end();
	   els_it != els_end; ++els_it, ++i) {
	 
	 double el_eta = els_it->Eta();
	 double el_phi = els_it->Phi();
	 double dR = deltaR(mu_eta, mu_phi, el_eta, el_phi);

	 if(dR < minDR) {
	   minDR = dR;
	   index = i;
	 }
       }

       // fill vector
       vector_mus_closestEle->push_back(index);
       vector_mus_eledr     ->push_back(minDR);
     }

     // store vectors
     iEvent.put(vector_mus_closestEle, "musclosestEle");
     iEvent.put(vector_mus_eledr     , "museledr"     );
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuToElsAssMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuToElsAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuToElsAssMaker);
