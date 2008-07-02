// -*- C++ -*-
//
// Package:    TrackToMuonAssMaker
// Class:      TrackToMuonAssMaker
// 
/**\class TrackToMuonAssMaker TrackToMuonAssMaker.cc CMS2/TrackToMuonAssMaker/src/TrackToMuonAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackToMuonAssMaker.cc,v 1.1 2008/07/02 03:32:45 jmuelmen Exp $
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
#include "CMS2/NtupleMaker/interface/TrackToMuonAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"

typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

TrackToMuonAssMaker::TrackToMuonAssMaker(const edm::ParameterSet& iConfig)
{
     produces<vector<int>   >("trkmusidx").setBranchAlias("trk_musidx");	// track index matched to muon
     produces<vector<float> >("trkmusdr" ).setBranchAlias("trk_musdr" );
}

void TrackToMuonAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>   > vector_trk_musidx(new vector<int>  );
     std::auto_ptr<vector<float> > vector_trk_musdr (new vector<float>);

     // get muon track p4's
     Handle<vector<LorentzVector> > mus_trk_p4_h;
     iEvent.getByLabel("muonMaker", "mustrkp4", mus_trk_p4_h);  

     // get track p4's
     Handle<vector<LorentzVector> > trks_p4_h;
     iEvent.getByLabel("trackMaker", "trkstrkp4", trks_p4_h);  

     for (vector<LorentzVector>::const_iterator track = trks_p4_h->begin(),
	    trks_end = trks_p4_h->end();
	  track != trks_end; ++track) { 

       // find track with min dR
       int trkidx = -999;
       double min_dR = 999;
       int i_track = 0;

       for (vector<LorentzVector>::const_iterator muon = mus_trk_p4_h->begin(),
	    muons_end = mus_trk_p4_h->end();
	    muon != muons_end; ++muon, ++i_track) {

	 const double deltaR = ROOT::Math::VectorUtil::DeltaR(*track, *muon);

	 if (deltaR < min_dR) {
	   min_dR = deltaR;
	   trkidx = i_track;
	 }
       }

       // fill vector
       vector_trk_musidx->push_back(trkidx);
       vector_trk_musdr ->push_back(min_dR );
     }

     // store vectors
     iEvent.put(vector_trk_musidx, "trkmusidx");
     iEvent.put(vector_trk_musdr , "trkmusdr" );
}

// ------------ method called once each job just before starting event loop  ------------
void 
TrackToMuonAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackToMuonAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackToMuonAssMaker);
