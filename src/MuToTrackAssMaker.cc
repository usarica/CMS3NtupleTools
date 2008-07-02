// -*- C++ -*-
//
// Package:    MuToTrackAssMaker
// Class:      MuToTrackAssMaker
// 
/**\class MuToTrackAssMaker MuToTrackAssMaker.cc CMS2/MuToTrackAssMaker/src/MuToTrackAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuToTrackAssMaker.cc,v 1.1 2008/07/02 02:20:18 jmuelmen Exp $
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
#include "CMS2/NtupleMaker/interface/MuToTrackAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"

typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

MuToTrackAssMaker::MuToTrackAssMaker(const edm::ParameterSet& iConfig)
{
     produces<vector<int>   >("mustrkidx").setBranchAlias("mus_trkidx");	// track index matched to muon
     produces<vector<float> >("mustrkdr" ).setBranchAlias("mus_trkdr" );
}

void MuToTrackAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>   > vector_mus_trkidx(new vector<int>  );
     std::auto_ptr<vector<float> > vector_mus_trkdr (new vector<float>);

     // get muon track p4's
     Handle<vector<LorentzVector> > mus_trk_p4_h;
     iEvent.getByLabel("muonMaker", "mustrkp4", mus_trk_p4_h);  

     // get track p4's
     Handle<vector<LorentzVector> > trks_p4_h;
     iEvent.getByLabel("trackMaker", "trkstrkp4", trks_p4_h);  

     vector<LorentzVector>::const_iterator muons_end = mus_trk_p4_h->end();

     for (vector<LorentzVector>::const_iterator muon = mus_trk_p4_h->begin(); 
	  muon != muons_end; ++muon) {

	  // find track with min dR
	  int trkidx = -999;
	  double min_dR = 999;
	  int i_track = 0;

	  for (vector<LorentzVector>::const_iterator track = trks_p4_h->begin(); 
	       track != trks_p4_h->end(); ++track, ++i_track) { 

	       const double deltaR = ROOT::Math::VectorUtil::DeltaR(*muon, *track);

	       if (deltaR < min_dR) {
		    min_dR = deltaR;
		    trkidx = i_track;
	       }
	  }

	  // fill vector
	  vector_mus_trkidx->push_back(trkidx);
	  vector_mus_trkdr ->push_back(min_dR);
     }

     // store vectors
     iEvent.put(vector_mus_trkidx, "mustrkidx");
     iEvent.put(vector_mus_trkdr , "mustrkdr" );
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuToTrackAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuToTrackAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuToTrackAssMaker);
