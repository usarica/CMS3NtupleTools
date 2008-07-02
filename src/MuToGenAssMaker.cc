// -*- C++ -*-
//
// Package:    MuToGenAssMaker
// Class:      MuToGenAssMaker
// 
/**\class MuToGenAssMaker MuToGenAssMaker.cc CMS2/NtupleMaker/src/MuToGenAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuToGenAssMaker.cc,v 1.1 2008/07/02 02:28:07 jmuelmen Exp $
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
#include "CMS2/NtupleMaker/interface/MuToGenAssMaker.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"


typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

MuToGenAssMaker::MuToGenAssMaker(const edm::ParameterSet& iConfig)
{
     produces<vector<int>           >("musmcid"      ).setBranchAlias("mus_mcid"      ); // muon matched to gen particle
     produces<vector<int>           >("musmcmotherid").setBranchAlias("mus_mcmotherid");
     produces<vector<int>           >("musmcidx"     ).setBranchAlias("mus_mcidx"     );
     produces<vector<LorentzVector> >("musmcp4"      ).setBranchAlias("mus_mcp4"      );

     //muonsInputTag        = iConfig.getParameter<edm::InputTag>("muonsInputTag"       );
     //genParticlesInputTag = iConfig.getParameter<edm::InputTag>("genParticlesInputTag");
}

void MuToGenAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>           > vector_mus_mcid      (new vector<int>          );
     std::auto_ptr<vector<int>           > vector_mus_mcmotherid(new vector<int>          );
     std::auto_ptr<vector<int>           > vector_mus_mcidx     (new vector<int>          );
     std::auto_ptr<vector<LorentzVector> > vector_mus_mcp4      (new vector<LorentzVector>);

     // get muons
     Handle<edm::View<reco::Muon> > muonhandle;
     iEvent.getByLabel("allLayer1TopMuons", muonhandle);      // change this in the future

     // get MC particle collection
     edm::Handle<reco::GenParticleCollection> genParticlesHandle;
     iEvent.getByLabel("genParticles", genParticlesHandle);

     for (edm::View<reco::Muon>::const_iterator muon = muonhandle->begin(),
	  muons_end = muonhandle->end();
	  muon != muons_end; ++muon) {

       //MC matching stuff
       int mcid = -999, mom_mcid = -999, genidx = -999;

       const reco::GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*muon, genParticlesHandle.product(), genidx);

       LorentzVector mc_p4(0,0,0,0);

       if(matchedGenParticle != 0) {
	 mcid                = matchedGenParticle->pdgId();
	 LorentzVector mc_p4 = matchedGenParticle->p4();
	 mom_mcid            = MCUtilities::motherID(*matchedGenParticle)->pdgId();
       }

       // fill vector
       vector_mus_mcid      ->push_back(mcid    );
       vector_mus_mcmotherid->push_back(mom_mcid);
       vector_mus_mcidx     ->push_back(genidx  );
       vector_mus_mcp4      ->push_back(mc_p4   );
     }

     // store vectors
     iEvent.put(vector_mus_mcid      , "musmcid"      );
     iEvent.put(vector_mus_mcmotherid, "musmcmotherid");
     iEvent.put(vector_mus_mcidx     , "musmcidx"     );
     iEvent.put(vector_mus_mcp4      , "musmcp4"      );
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuToGenAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuToGenAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuToGenAssMaker);
