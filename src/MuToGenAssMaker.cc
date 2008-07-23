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
// $Id: MuToGenAssMaker.cc,v 1.4 2008/07/23 05:35:01 fgolf Exp $
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
using namespace std;
using namespace edm;

MuToGenAssMaker::MuToGenAssMaker(const edm::ParameterSet& iConfig)
{
     produces<vector<int>           >("musmcid"      ).setBranchAlias("mus_mc_id"      ); // muon matched to gen particle
     produces<vector<int>           >("musmcmotherid").setBranchAlias("mus_mc_motherid");
     produces<vector<int>           >("musmcidx"     ).setBranchAlias("mus_mcidx"      );
     produces<vector<LorentzVector> >("musmcp4"      ).setBranchAlias("mus_mcp4"       );
     produces<vector<double>        >("musmcdr"      ).setBranchAlias("mus_mcdr"       );
}

void MuToGenAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>           > vector_mus_mc_id      (new vector<int>          );
     std::auto_ptr<vector<int>           > vector_mus_mc_motherid(new vector<int>          );
     std::auto_ptr<vector<int>           > vector_mus_mcidx      (new vector<int>          );
     std::auto_ptr<vector<LorentzVector> > vector_mus_mcp4       (new vector<LorentzVector>);
     std::auto_ptr<vector<double>        > vector_mus_mcdr       (new vector<double>       );

     // get muons
     Handle<vector<LorentzVector> > mus_p4_h;
     iEvent.getByLabel("muonMaker", "musp4", mus_p4_h);

     // get MC particles
     Handle<vector<LorentzVector> > genps_p4_h;
     iEvent.getByLabel("genMaker", "genpsp4", genps_p4_h);

     Handle<vector<int> > genps_id_h;
     iEvent.getByLabel("genMaker", "genpsid", genps_id_h);

     Handle<vector<int> > genps_id_mother_h;
     iEvent.getByLabel("genMaker", "genpsidmother", genps_id_mother_h);

     for (vector<LorentzVector>::const_iterator muon = mus_p4_h->begin(),
	  mus_end = mus_p4_h->end();
	  muon != mus_end; ++muon) { 

       //MC matching stuff
       int mcid = -999, mom_mcid = -999, genidx = -999;
       LorentzVector mc_p4(0,0,0,0);

       int genp = 0; // gen particle counter
       double min_dR = 999;

       for (vector<LorentzVector>::const_iterator genps = genps_p4_h->begin(),
	    genps_end = genps_p4_h->end();
	    genps != genps_end; ++genps, ++genp) { 

	 const double deltaR = ROOT::Math::VectorUtil::DeltaR(*muon, *genps);

	 if (deltaR < min_dR) {
	   min_dR   = deltaR;
	   mcid     = (*genps_id_h)[genp];
	   genidx   = genp;
	   mom_mcid = (*genps_id_mother_h)[genp];
	   mc_p4    = (*genps_p4_h)[genp];

	 }	 

       }

       // fill vector
       vector_mus_mc_id      ->push_back(mcid    );
       vector_mus_mc_motherid->push_back(mom_mcid);
       vector_mus_mcidx      ->push_back(genidx  );
       vector_mus_mcp4       ->push_back(mc_p4   );
       vector_mus_mcdr       ->push_back(min_dR  );
     }

     // store vectors
     iEvent.put(vector_mus_mc_id      , "musmcid"      );
     iEvent.put(vector_mus_mc_motherid, "musmcmotherid");
     iEvent.put(vector_mus_mcidx      , "musmcidx"     );
     iEvent.put(vector_mus_mcp4       , "musmcp4"      );
     iEvent.put(vector_mus_mcdr       , "musmcdr"      );
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
