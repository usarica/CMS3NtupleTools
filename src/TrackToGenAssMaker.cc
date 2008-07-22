// -*- C++ -*-
//
// Package:    TrackToGenAssMaker
// Class:      TrackToGenAssMaker
// 
/**\class TrackToGenAssMaker TrackToGenAssMaker.cc CMS2/NtupleMaker/src/TrackToGenAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackToGenAssMaker.cc,v 1.3 2008/07/22 20:12:20 fgolf Exp $
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
#include "CMS2/NtupleMaker/interface/TrackToGenAssMaker.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

typedef math::XYZTLorentzVector LorentzVector;
using std::vector;
using namespace edm;

TrackToGenAssMaker::TrackToGenAssMaker(const edm::ParameterSet& iConfig)
{
     produces<vector<int>           >("trkmcid"      ).setBranchAlias("trk_mc_id"      ); // track matched to gen particle
     produces<vector<int>           >("trkmcmotherid").setBranchAlias("trk_mc_motherid");
     produces<vector<int>           >("trkmcidx"     ).setBranchAlias("trk_mcidx"      );
     produces<vector<LorentzVector> >("trkmcp4"      ).setBranchAlias("trk_mcp4"       );
     produces<vector<double>        >("trkmcdr"      ).setBranchAlias("trk_mcdr"       );
}

void TrackToGenAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>           > vector_trk_mc_id      (new vector<int>          );
     std::auto_ptr<vector<int>           > vector_trk_mc_motherid(new vector<int>          );
     std::auto_ptr<vector<int>           > vector_trk_mcidx      (new vector<int>          );
     std::auto_ptr<vector<LorentzVector> > vector_trk_mcp4       (new vector<LorentzVector>);
     std::auto_ptr<vector<double>        > vector_trk_mcdr       (new vector<double>       );

     // get tracks
     Handle<vector<LorentzVector> > trks_p4_h;
     iEvent.getByLabel("trackMaker", "trkstrkp4", trks_p4_h);

     // get MC particles
     Handle<vector<LorentzVector> > genps_p4_h;
     iEvent.getByLabel("genMaker", "genpsp4", genps_p4_h);

     Handle<vector<int> > genps_id_h;
     iEvent.getByLabel("genMaker", "genpsid", genps_id_h);

     Handle<vector<int> > genps_id_mother_h;
     iEvent.getByLabel("genMaker", "genpsidmother", genps_id_mother_h);

     for (vector<LorentzVector>::const_iterator track = trks_p4_h->begin(),
	  trks_end = trks_p4_h->end();
	  track != trks_end; ++track) { 

       //MC matching stuff
       int mcid = -999, mom_mcid = -999, genidx = -999;
       LorentzVector mc_p4(0,0,0,0);

       int genp = 0; // gen particle counter
       double min_dR = 999;

       for (vector<LorentzVector>::const_iterator genps = genps_p4_h->begin(),
	    genps_end = genps_p4_h->end();
	    genps != genps_end; ++genp) { 

	 const double deltaR = ROOT::Math::VectorUtil::DeltaR(*track, *genps);

	 if (deltaR < min_dR) {
	   min_dR   = deltaR;
	   mcid     = (*genps_id_h)[genp];
	   genidx   = genp;
	   mom_mcid = (*genps_id_mother_h)[genp];
	   mc_p4    = (*genps_p4_h)[genp];
	 }	 
       }

       // fill vector
       vector_trk_mc_id      ->push_back(mcid    );
       vector_trk_mc_motherid->push_back(mom_mcid);
       vector_trk_mcidx      ->push_back(genidx  );
       vector_trk_mcp4       ->push_back(mc_p4   );
       vector_trk_mcdr       ->push_back(min_dR  );
     }

     // store vectors
     iEvent.put(vector_trk_mc_id      , "trkmcid"      );
     iEvent.put(vector_trk_mc_motherid, "trkmcmotherid");
     iEvent.put(vector_trk_mcidx      , "trkmcidx"     );
     iEvent.put(vector_trk_mcp4       , "trkmcp4"      );
     iEvent.put(vector_trk_mcdr       , "trkmcdr"      );
}

// ------------ method called once each job just before starting event loop  ------------
void 
TrackToGenAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackToGenAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackToGenAssMaker);
