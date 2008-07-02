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
// $Id: TrackToGenAssMaker.cc,v 1.1 2008/07/02 03:32:45 jmuelmen Exp $
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

TrackToGenAssMaker::TrackToGenAssMaker(const edm::ParameterSet& iConfig)
{
     produces<vector<int>           >("trkmcid"      ).setBranchAlias("trk_mcid"      ); // track matched to gen particle
     produces<vector<int>           >("trkmcmotherid").setBranchAlias("trk_mcmotherid");
     produces<vector<int>           >("trkmcidx"     ).setBranchAlias("trk_mcidx"     );
     produces<vector<LorentzVector> >("trkmcp4"      ).setBranchAlias("trk_mcp4"      );

     //tracksInputTag       = iConfig.getParameter<InputTag>("tracksInputTag"      );
     //genParticlesInputTag = iConfig.getParameter<InputTag>("genParticlesInputTag");
}

void TrackToGenAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>           > vector_trk_mcid      (new vector<int>          );
     std::auto_ptr<vector<int>           > vector_trk_mcmotherid(new vector<int>          );
     std::auto_ptr<vector<int>           > vector_trk_mcidx     (new vector<int>          );
     std::auto_ptr<vector<LorentzVector> > vector_trk_mcp4      (new vector<LorentzVector>);

     // get tracks
     Handle<edm::View<reco::Track> > trackhandle;
     iEvent.getByLabel("ctfWithMaterialTracks", trackhandle);

     // get MC particle collection
     edm::Handle<reco::GenParticleCollection> genParticlesHandle;
     iEvent.getByLabel("genParticles", genParticlesHandle);

     for (edm::View<reco::Track>::const_iterator track = trackhandle->begin(),
	  tracks_end = trackhandle->end();
	  track != tracks_end; ++track) {

       //MC matching stuff

       int mcid = -999, mom_mcid = -999, genidx = -999;

       const reco::GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*track, genParticlesHandle.product(), genidx);

       LorentzVector mc_p4(0,0,0,0);

       if(matchedGenParticle != 0) {
	 mcid                = matchedGenParticle->pdgId();
	 LorentzVector mc_p4 = matchedGenParticle->p4();
	 mom_mcid            = MCUtilities::motherID(*matchedGenParticle)->pdgId();
       }

       // fill vector
       vector_trk_mcid      ->push_back(mcid    );
       vector_trk_mcmotherid->push_back(mom_mcid);
       vector_trk_mcidx     ->push_back(genidx  );
       vector_trk_mcp4      ->push_back(mc_p4   );
     }

     // store vectors
     iEvent.put(vector_trk_mcid      , "trkmcid"      );
     iEvent.put(vector_trk_mcmotherid, "trkmcmotherid");
     iEvent.put(vector_trk_mcidx     , "trkmcidx"     );
     iEvent.put(vector_trk_mcp4      , "trkmcp4"      );
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
