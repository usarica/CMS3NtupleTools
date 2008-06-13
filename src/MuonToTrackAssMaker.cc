// -*- C++ -*-
//
// Package:    MuonToTrackAssMaker
// Class:      MuonToTrackAssMaker
// 
/**\class MuonToTrackAssMaker MuonToTrackAssMaker.cc CMS2/MuonToTrackAssMaker/src/MuonToTrackAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonToTrackAssMaker.cc,v 1.2 2008/06/13 02:26:22 jmuelmen Exp $
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
#include "CMS2/NtupleMaker/interface/MuonToTrackAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"

typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

MuonToTrackAssMaker::MuonToTrackAssMaker(const edm::ParameterSet& iConfig)
{
     produces<vector<int> >("mustrkidx").setBranchAlias("mus_trkidx");	// track index matched to muon
}

void MuonToTrackAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int> > vector_mus_trkidx(new vector<int>);        
     // get muon p4's
     Handle<vector<LorentzVector> > mus_p4_h;
     iEvent.getByLabel("muonMaker", "musp4", mus_p4_h);  
     // get track p4's
     Handle<vector<LorentzVector> > trks_p4_h;
     iEvent.getByLabel("trackMaker", "trkstrkp4", trks_p4_h);  
     vector<LorentzVector>::const_iterator muons_end = mus_p4_h->end();
     for (vector<LorentzVector>::const_iterator muon = mus_p4_h->begin(); 
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
     }
     // store vectors
     iEvent.put(vector_mus_trkidx, "mustrkidx");
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonToTrackAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonToTrackAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonToTrackAssMaker);
