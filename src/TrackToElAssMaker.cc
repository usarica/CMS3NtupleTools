// -*- C++ -*-
//
// Package:    TrackToElAssMaker
// Class:      TrackToElAssMaker
// 
/**\class TrackToElAssMaker TrackToElAssMaker.cc CMS2/NtupleMaker/src/TrackToElAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackToElAssMaker.cc,v 1.1 2008/07/02 03:32:45 jmuelmen Exp $
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
#include "CMS2/NtupleMaker/interface/TrackToElAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

TrackToElAssMaker::TrackToElAssMaker(const edm::ParameterSet& iConfig)
{
     produces<vector<int>   >("trkelsidx").setBranchAlias("trk_elsidx");	// track index matched to electron
     produces<vector<float> >("trkelsdr" ).setBranchAlias("trk_elsdr" );
}

void TrackToElAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  // make vectors to hold the information
  std::auto_ptr<vector<int>   > vector_trk_elsidx(new vector<int>  );
  std::auto_ptr<vector<float> > vector_trk_elsdr (new vector<float>);

  // get electron p4's
  Handle<vector<LorentzVector> > els_p4_h;
  iEvent.getByLabel("electronMaker", "elsp4", els_p4_h);  

  // get track p4's
  Handle<vector<LorentzVector> > trks_p4_h;
  iEvent.getByLabel("trackMaker", "trkstrkp4", trks_p4_h);

  // loop over tracks
  for(vector<LorentzVector>::const_iterator trk_it = trks_p4_h->begin(),
	trk_end = trks_p4_h->end();
      trk_it != trk_end; ++trk_it) {
	 
    double trk_eta = trk_it->Eta();
    double trk_phi = trk_it->Phi();
       
    double minDR   = 999;
    unsigned int i = 0;
    int index      = -999;
     
    //loop over electrons and find the closest track
    for(vector<LorentzVector>::const_iterator els_it = els_p4_h->begin(),
	  els_end = els_p4_h->end();
	els_it != els_end; ++els_it, ++i) {
       
      double el_eta = els_it->Eta();
      double el_phi = els_it->Phi();

      double dR = deltaR(trk_eta, trk_phi, el_eta, el_phi);

      if(dR < minDR) {
	minDR = dR;
	index = i;
      }
    }

    // fill vector
    vector_trk_elsidx->push_back(index);
    vector_trk_elsdr ->push_back(minDR);
  }

  // store vectors
  iEvent.put(vector_trk_elsidx, "trkelsidx");
  iEvent.put(vector_trk_elsdr , "trkelsdr" );
}

// ------------ method called once each job just before starting event loop  ------------
void 
TrackToElAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackToElAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackToElAssMaker);
