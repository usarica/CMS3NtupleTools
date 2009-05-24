// -*- C++ -*-
//
// Package:    PATMuonMaker
// Class:      PATMuonMaker
// 
/**\class PATMuonMaker PATMuonMaker.cc CMS2/NtupleMaker/src/PATMuonMaker.cc

Description: copy additional PAT muon variables in simple data structures into the EDM event tree

 Implementation:
     - take PAT muons
     - extract and fill variables
*/
//
// Original Author:  Frank Golf
// Thu Jun 25 16:39:55 UTC 2008
// $Id: PATMuonMaker.cc,v 1.8 2009/05/24 19:36:51 kalavase Exp $
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/PATMuonMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Flags.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;

PATMuonMaker::PATMuonMaker(const edm::ParameterSet& iConfig) {

  // product of this EDProducer
  produces<vector<LorentzVector>  >   ("muspatp4"          ).setBranchAlias("mus_pat_p4"          );
  produces<vector<float>          >   ("muspattrackIso"    ).setBranchAlias("mus_pat_trackIso"    );
  produces<vector<float>          >   ("muspatcaloIso"     ).setBranchAlias("mus_pat_caloIso"     );
  produces<vector<float>          >   ("muspatecalIso"     ).setBranchAlias("mus_pat_ecalIso"     );
  produces<vector<float>          >   ("muspathcalIso"     ).setBranchAlias("mus_pat_hcalIso"     );
  //deposit in the veto cone - gotten from the IsoDeposits
  produces<vector<float>          >   ("muspatvetoDep"     ).setBranchAlias("mus_pat_vetoDep"     );
  produces<vector<float>          >   ("muspatecalvetoDep" ).setBranchAlias("mus_pat_ecalvetoDep" );
  produces<vector<float>          >   ("muspathcalvetoDep" ).setBranchAlias("mus_pat_hcalvetoDep" );
  produces<vector<int>            >   ("muspatgenID"       ).setBranchAlias("mus_pat_genID"       );
  produces<vector<int>            >   ("muspatgenMotherID" ).setBranchAlias("mus_pat_genMotherID" );
  produces<vector<uint32_t>       >   ("muspatflag"        ).setBranchAlias("mus_pat_flag"        );
  produces<vector<LorentzVector>  >   ("muspatgenP4"       ).setBranchAlias("mus_pat_genP4"       );
  produces<vector<LorentzVector>  >   ("muspatgenMotherP4" ).setBranchAlias("mus_pat_genMotherP4" );

  // parameters from configuration
  patMuonsInputTag_  = iConfig.getParameter<edm::InputTag>("patMuonsInputTag"  );
  recoMuonsInputTag_ = iConfig.getParameter<edm::InputTag>("recoMuonsInputTag" );
}

PATMuonMaker::~PATMuonMaker() {}

// ------------ method called to produce the data  ------------
void PATMuonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // get jet collection
  edm::Handle<vector<pat::Muon> > patMuonHandle;
  iEvent.getByLabel(patMuonsInputTag_, patMuonHandle);
  vector<pat::Muon> v_patMuons = *(patMuonHandle.product());

  edm::Handle<vector<reco::Muon> > recoMuonHandle;
  iEvent.getByLabel(recoMuonsInputTag_, recoMuonHandle);
  vector<reco::Muon> v_recoMuons = *(recoMuonHandle.product());
  
  //make sure that the PAT and recoMuon collections are indeed aligned 
  MatchUtilities::alignRecoPatMuonCollections(v_recoMuons, v_patMuons);

  // create containers
  auto_ptr<vector<LorentzVector>     >   mus_pat_p4          ( new vector<LorentzVector>     );
  auto_ptr<vector<float>             >   mus_pat_trackIso    ( new vector<float>             );
  auto_ptr<vector<float>             >   mus_pat_caloIso     ( new vector<float>             );
  auto_ptr<vector<float>             >   mus_pat_ecalIso     ( new vector<float>             );
  auto_ptr<vector<float>             >   mus_pat_hcalIso     ( new vector<float>             );
  auto_ptr<vector<float>             >   mus_pat_vetoDep     ( new vector<float>             );
  auto_ptr<vector<float>             >   mus_pat_ecalvetoDep ( new vector<float>             );
  auto_ptr<vector<float>             >   mus_pat_hcalvetoDep ( new vector<float>             );
  auto_ptr<vector<int>               >   mus_pat_genID       ( new vector<int>               );
  auto_ptr<vector<int>               >   mus_pat_genMotherID ( new vector<int>               ); 
  auto_ptr<vector<uint32_t>          >   mus_pat_flag        ( new vector<uint32_t>          ); 
  auto_ptr<vector<LorentzVector>     >   mus_pat_genP4       ( new vector<LorentzVector>     );
  auto_ptr<vector<LorentzVector>     >   mus_pat_genMotherP4 ( new vector<LorentzVector>     );
  
  // loop over top muons and fill containers
  for ( std::vector<pat::Muon>::const_iterator patmu_it = v_patMuons.begin(),
	patmu_end = v_patMuons.end();
	patmu_it != patmu_end; patmu_it++) {
    
    GenParticle gen(patmu_it->genLepton() ? *patmu_it->genLepton() : 
		    reco::GenParticle(0, reco::Particle::LorentzVector(0, 0, 0, 0), reco::Particle::Point(0,0,0), 0, 0, true));
    const GenParticle *gen_mom = MCUtilities::motherID(gen);


    //the veto cone size is 0.07 in the ecal and 0.1 in the hcal 
    const reco::IsoDeposit *ecalIsoDep = patmu_it->ecalIsoDeposit();
    const reco::IsoDeposit *hcalIsoDep = patmu_it->hcalIsoDeposit();
    
    mus_pat_p4          ->push_back( patmu_it->p4()             );
    mus_pat_trackIso    ->push_back( patmu_it->trackIso()       );
    mus_pat_vetoDep     ->push_back( ecalIsoDep->candEnergy()
				     + hcalIsoDep->candEnergy() );
    mus_pat_ecalvetoDep ->push_back( ecalIsoDep->candEnergy()   );
    mus_pat_hcalvetoDep ->push_back( hcalIsoDep->candEnergy()   );
    mus_pat_ecalIso     ->push_back( patmu_it->ecalIso()        );
    mus_pat_hcalIso     ->push_back( patmu_it->hcalIso()        );
    mus_pat_caloIso     ->push_back( patmu_it->caloIso()        );
    mus_pat_genID       ->push_back( gen.pdgId()            );
    mus_pat_genMotherID ->push_back( gen_mom->pdgId()       );
    mus_pat_flag        ->push_back( patmu_it->status()     );
    mus_pat_genP4       ->push_back( gen.p4()               );
    mus_pat_genMotherP4 ->push_back( gen_mom->p4()          );
    

  }
  

  // put containers into event
  iEvent.put(mus_pat_p4,           "muspatp4"              );
  iEvent.put(mus_pat_trackIso,     "muspattrackIso"        );
  iEvent.put(mus_pat_vetoDep,      "muspatvetoDep"         );
  iEvent.put(mus_pat_ecalvetoDep,  "muspatecalvetoDep"     );
  iEvent.put(mus_pat_hcalvetoDep,  "muspathcalvetoDep"     );
  iEvent.put(mus_pat_caloIso,      "muspatcaloIso"         );
  iEvent.put(mus_pat_ecalIso,      "muspatecalIso"         );
  iEvent.put(mus_pat_hcalIso,      "muspathcalIso"         );
  
  iEvent.put(mus_pat_genID,        "muspatgenID"           );
  iEvent.put(mus_pat_genMotherID,  "muspatgenMotherID"     );
  iEvent.put(mus_pat_flag,         "muspatflag"            );
  iEvent.put(mus_pat_genP4,        "muspatgenP4"           );
  iEvent.put(mus_pat_genMotherP4,  "muspatgenMotherP4"     );
}

// ------------ method called once each job just before starting event loop  ------------
void 
PATMuonMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATMuonMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATMuonMaker);
