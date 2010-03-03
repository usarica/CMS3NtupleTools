// -*- C++ -*-
//
// Package:    TrackToMuonAssMaker
// Class:      TrackToMuonAssMaker
// 
/**\class TrackToMuonAssMaker TrackToMuonAssMaker.cc CMS2/TrackToMuonAssMaker/src/TrackToMuonAssMaker.cc

   Description: <one line class summary>

   Implementation:
   Code to associate a track to a muon. Since a GlobalMuon or TrackerMuon
   will always have a Si track, we want to keep the dR requirement very 
   tight, so that we always associate the track to a muon only if it 
   actually belongs to a muon
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackToMuonAssMaker.cc,v 1.8 2010/03/03 04:24:06 kalavase Exp $
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

typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

TrackToMuonAssMaker::TrackToMuonAssMaker(const edm::ParameterSet& iConfig) {


  produces<vector<int>   >("trkmusidx").setBranchAlias("trk_musidx");	// track index matched to muon
     
  m_minDR_      = iConfig.getParameter<double>("minDR");
  musInputTag_  = iConfig.getParameter<edm::InputTag>("musInputTag");
  trksInputTag_ = iConfig.getParameter<edm::InputTag>("trksInputTag");
}

void TrackToMuonAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  // make vectors to hold the information
  std::auto_ptr<vector<int>   > vector_trk_musidx(new vector<int>  );

  // get muon track p4's
  Handle<vector<LorentzVector> > mus_trk_p4_h;
  iEvent.getByLabel(musInputTag_.label(), "mustrkp4", mus_trk_p4_h);

  if( !mus_trk_p4_h.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve mus_trk_p4";
    edm::LogInfo("OutputInfo") << " TrackToMuonAssMaker cannot continue...!";
    return;
  }

  //get the muon type
  Handle<vector<int> > mus_type_h;
  iEvent.getByLabel(musInputTag_.label(), "mustype", mus_type_h);

  if( !mus_type_h.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve mus_type";
    edm::LogInfo("OutputInfo") << " TrackToMuonAssMaker cannot continue...!";
    return;
  }

  const vector<int> *mus_type = mus_type_h.product();
     
  // get track p4's
  Handle<vector<LorentzVector> > trks_p4_h;
  iEvent.getByLabel(trksInputTag_.label(), "trkstrkp4", trks_p4_h);  

  if( !trks_p4_h.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve trks_p4";
    edm::LogInfo("OutputInfo") << " TrackToMuonAssMaker cannot continue...!";
    return;
  }  
     
  for (vector<LorentzVector>::const_iterator track = trks_p4_h->begin(), trks_end = trks_p4_h->end(); track != trks_end; ++track) { 
	  
    double minDR = m_minDR_;
    int trkidx   = -9999;
    int i_track  = 0;

    for (vector<LorentzVector>::const_iterator muon = mus_trk_p4_h->begin(), muons_end = mus_trk_p4_h->end(); muon != muons_end; ++muon, ++i_track) {
	 
      if( !( mus_type->at(i_track) & 1<<1 ) && !( mus_type->at(i_track) & 1<<2 ) ) continue;       //do the association ONLY if it is either a tracker muon or global muon
	 
      const double deltaR = ROOT::Math::VectorUtil::DeltaR(*track, *muon);
	 
      if (deltaR < minDR) {
	minDR  = deltaR;
	trkidx = i_track;
      }
    }

    // fill vector
    vector_trk_musidx->push_back( trkidx );
  }

  // store vectors
  iEvent.put(vector_trk_musidx, "trkmusidx" );
}

// ------------ method called once each job just before starting event loop  ------------
void TrackToMuonAssMaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackToMuonAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackToMuonAssMaker);
