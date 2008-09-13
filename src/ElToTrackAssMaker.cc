// -*- C++ -*-
//
// Package:    ElToTrackAssMaker
// Class:      ElToTrackAssMaker
// 
/**\class ElToTrackAssMaker ElToTrackAssMaker.cc CMS2/NtupleMaker/src/ElToTrackAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElToTrackAssMaker.cc,v 1.2 2008/09/13 08:07:23 jmuelmen Exp $
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
#include "CMS2/NtupleMaker/interface/ElToTrackAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

ElToTrackAssMaker::ElToTrackAssMaker(const edm::ParameterSet& iConfig)
     : m_minDR(iConfig.getParameter<double>("minDR"))
{
     produces<vector<int>   >("elstrkidx").setBranchAlias("els_trkidx");	// track index matched to electron
     produces<vector<float> >("elstrkdr" ).setBranchAlias("els_trkdr" );
}

void ElToTrackAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>   > vector_els_trkidx(new vector<int>  );
     std::auto_ptr<vector<float> > vector_els_trkdr (new vector<float>);

     // get electron p4's
     Handle<vector<LorentzVector> > els_p4_h;
     iEvent.getByLabel("electronMaker", "elsp4", els_p4_h);  

     // get track p4's
     Handle<vector<LorentzVector> > trks_p4_h;
     iEvent.getByLabel("trackMaker", "trkstrkp4", trks_p4_h);
     
     //loop over electrons and find the closest track
     for(vector<LorentzVector>::const_iterator els_it = els_p4_h->begin(),
	 els_end = els_p4_h->end();
	 els_it != els_end; els_it++) {
       
       double el_eta = els_it->Eta();
       double el_phi = els_it->Phi();
       
       double minDR = m_minDR;
       unsigned int i = 0;
       int index      = -999;

       // loop over tracks
       for(vector<LorentzVector>::const_iterator trk_it = trks_p4_h->begin(),
	   trk_end = trks_p4_h->end();
	   trk_it != trk_end; ++trk_it, ++i) {
	 
	 double trk_eta = trk_it->Eta();
	 double trk_phi = trk_it->Phi();

	 double dR = deltaR(el_eta, el_phi, trk_eta, trk_phi);

	 if(dR < minDR) {
	   minDR = dR;
	   index = i;
	 }
       }

       // fill vector
       vector_els_trkidx->push_back(index);
       vector_els_trkdr ->push_back(minDR);
     }

     // store vectors
     iEvent.put(vector_els_trkidx, "elstrkidx");
     iEvent.put(vector_els_trkdr , "elstrkdr" );
}

// ------------ method called once each job just before starting event loop  ------------
void 
ElToTrackAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElToTrackAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElToTrackAssMaker);
