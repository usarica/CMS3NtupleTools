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
// $Id: PATMuonMaker.cc,v 1.14 2010/03/18 02:13:17 kalavase Exp $
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

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;

PATMuonMaker::PATMuonMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // product of this EDProducer
  produces<vector<LorentzVector>  >   (branchprefix+"p4"          ).setBranchAlias(aliasprefix_+"_p4"          );
  produces<vector<float>          >   (branchprefix+"trackIso"    ).setBranchAlias(aliasprefix_+"_trackIso"    );
  produces<vector<float>          >   (branchprefix+"caloIso"     ).setBranchAlias(aliasprefix_+"_caloIso"     );
  produces<vector<float>          >   (branchprefix+"ecalIso"     ).setBranchAlias(aliasprefix_+"_ecalIso"     );
  produces<vector<float>          >   (branchprefix+"hcalIso"     ).setBranchAlias(aliasprefix_+"_hcalIso"     );
  //deposit in the veto cone - gotten from the IsoDeposits
  produces<vector<float>          >   (branchprefix+"trckvetoDep" ).setBranchAlias(aliasprefix_+"_trckvetoDep" );
  produces<vector<float>          >   (branchprefix+"calovetoDep" ).setBranchAlias(aliasprefix_+"_calovetoDep" );
  produces<vector<float>          >   (branchprefix+"ecalvetoDep" ).setBranchAlias(aliasprefix_+"_ecalvetoDep" );
  produces<vector<float>          >   (branchprefix+"hcalvetoDep" ).setBranchAlias(aliasprefix_+"_hcalvetoDep" );
  produces<vector<int>            >   (branchprefix+"genID"       ).setBranchAlias(aliasprefix_+"_genID"       );
  produces<vector<int>            >   (branchprefix+"genMotherID" ).setBranchAlias(aliasprefix_+"_genMotherID" );
  produces<vector<uint32_t>       >   (branchprefix+"flag"        ).setBranchAlias(aliasprefix_+"_flag"        );
  produces<vector<LorentzVector>  >   (branchprefix+"genP4"       ).setBranchAlias(aliasprefix_+"_genP4"       );
  produces<vector<LorentzVector>  >   (branchprefix+"genMotherP4" ).setBranchAlias(aliasprefix_+"_genMotherP4" );

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
  auto_ptr<vector<float>             >   mus_pat_trckvetoDep ( new vector<float>             ); //energy in track veto cone
  auto_ptr<vector<float>             >   mus_pat_calovetoDep ( new vector<float>             ); //sum of ecal+hcal veto cones
  auto_ptr<vector<float>             >   mus_pat_ecalvetoDep ( new vector<float>             ); //energy in ecal veto cone
  auto_ptr<vector<float>             >   mus_pat_hcalvetoDep ( new vector<float>             ); //energy in hcal veto cone
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
    // and 0.01 in the tracker
    const reco::IsoDeposit *trckIsoDep = patmu_it->trackIsoDeposit();
    const reco::IsoDeposit *ecalIsoDep = patmu_it->ecalIsoDeposit();
    const reco::IsoDeposit *hcalIsoDep = patmu_it->hcalIsoDeposit();

    mus_pat_p4          ->push_back( LorentzVector( patmu_it->p4() ) );
    mus_pat_trackIso    ->push_back( patmu_it->userIsolation(pat::TrackIso) );
    mus_pat_caloIso     ->push_back( patmu_it->userIsolation(pat::EcalIso) + patmu_it->userIsolation(pat::HcalIso) );
    mus_pat_ecalIso     ->push_back( patmu_it->userIsolation(pat::EcalIso)  );
    mus_pat_hcalIso     ->push_back( patmu_it->userIsolation(pat::HcalIso)  );
    mus_pat_trckvetoDep ->push_back( trckIsoDep->candEnergy() );
    mus_pat_calovetoDep ->push_back( ecalIsoDep->candEnergy() + hcalIsoDep->candEnergy() );
    mus_pat_ecalvetoDep ->push_back( ecalIsoDep->candEnergy() );
    mus_pat_hcalvetoDep ->push_back( hcalIsoDep->candEnergy() );
    mus_pat_genID       ->push_back( gen.pdgId()              );
    mus_pat_genMotherID ->push_back( gen_mom->pdgId()         );
    mus_pat_flag        ->push_back( patmu_it->status()       );
    mus_pat_genP4       ->push_back( LorentzVector( gen.p4() )      );
    mus_pat_genMotherP4 ->push_back( LorentzVector( gen_mom->p4() ) );
    

  }
  

  // put containers into event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(mus_pat_p4,           branchprefix+"p4"              );
  iEvent.put(mus_pat_trackIso,     branchprefix+"trackIso"        );
  iEvent.put(mus_pat_caloIso,      branchprefix+"caloIso"         );
  iEvent.put(mus_pat_ecalIso,      branchprefix+"ecalIso"         );
  iEvent.put(mus_pat_hcalIso,      branchprefix+"hcalIso"         );
  iEvent.put(mus_pat_trckvetoDep,  branchprefix+"trckvetoDep"     );
  iEvent.put(mus_pat_calovetoDep,  branchprefix+"calovetoDep"     );
  iEvent.put(mus_pat_ecalvetoDep,  branchprefix+"ecalvetoDep"     );
  iEvent.put(mus_pat_hcalvetoDep,  branchprefix+"hcalvetoDep"     );
  
  iEvent.put(mus_pat_genID,        branchprefix+"genID"           );
  iEvent.put(mus_pat_genMotherID,  branchprefix+"genMotherID"     );
  iEvent.put(mus_pat_flag,         branchprefix+"flag"            );
  iEvent.put(mus_pat_genP4,        branchprefix+"genP4"           );
  iEvent.put(mus_pat_genMotherP4,  branchprefix+"genMotherP4"     );
}

// ------------ method called once each job just before starting event loop  ------------
void 
PATMuonMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATMuonMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATMuonMaker);
