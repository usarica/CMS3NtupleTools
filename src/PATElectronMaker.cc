// -*- C++ -*-
//
// Package:    PATElectronMaker
// Class:      PATElectronMaker
// 
/**\class PATElectronMaker PATElectronMaker.cc CMS2/NtupleMaker/src/PATElectronMaker.cc

Description: copy additional PAT jet variables in simple data structures into the EDM event tree

Implementation:
- take PAT jets
- extract and fill variables
*/
//
// Original Author:  Puneeth Kalavase
// Thu Jun 12 22:55:46 UTC 2008
// $Id: PATElectronMaker.cc,v 1.6 2009/05/24 19:36:13 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/PATElectronMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Flags.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;



PATElectronMaker::PATElectronMaker(const edm::ParameterSet& iConfig) {

  // product of this EDProducer
  produces<vector<LorentzVector>  >   ("elspatp4").setBranchAlias("els_pat_p4");
  produces<vector<int>            >   ("elspatgenID"               ).setBranchAlias("els_pat_genID"            );
  produces<vector<int>            >   ("elspatgenMotherID"         ).setBranchAlias("els_pat_genMotherID"      );
  produces<vector<uint32_t>       >   ("elspatflag"                ).setBranchAlias("els_pat_flag"             );
  produces<vector<float>          >   ("elspattrackIso"            ).setBranchAlias("els_pat_trackIso"         );
  produces<vector<float>          >   ("elspatcaloIso"             ).setBranchAlias("els_pat_caloIso"          );
  produces<vector<float>          >   ("elspatecalIso"             ).setBranchAlias("els_pat_ecalIso"          );
  produces<vector<float>          >   ("elspathcalIso"             ).setBranchAlias("els_pat_hcalIso"          );
  produces<vector<float>          >   ("elspatrobustLooseId"       ).setBranchAlias("els_pat_robustLooseId"    );
  produces<vector<float>          >   ("elspatrobustTightId"       ).setBranchAlias("els_pat_robustTightId"    );
  produces<vector<float>          >   ("elspatlooseId"             ).setBranchAlias("els_pat_looseId"          );
  produces<vector<float>          >   ("elspattightId"             ).setBranchAlias("els_pat_tightId"          );
  produces<vector<float>          >   ("elspatrobustHighEnergy"    ).setBranchAlias("els_pat_robustHighEnergy" );
  produces<vector<LorentzVector>  >   ("elspatgenP4"               ).setBranchAlias("els_pat_genP4"            );
  produces<vector<LorentzVector>  >   ("elspatgenMotherP4"         ).setBranchAlias("els_pat_genMotherP4"      );
  produces<vector<float>          >   ("elspatsigmaEtaEta"         ).setBranchAlias("els_pat_sigmaEtaEta"      );
  produces<vector<float>          >   ("elspatsigmaIEtaIEta"       ).setBranchAlias("els_pat_sigmaIEtaIEta"    );
  produces<vector<float>          >   ("elspatscE1x5"              ).setBranchAlias("els_pat_scE1x5"           );
  produces<vector<float>          >   ("elspatscE2x5Max"           ).setBranchAlias("els_pat_scE2x5Max"        );
  produces<vector<float>          >   ("elspatscE5x5"              ).setBranchAlias("els_pat_scE5x5"           );

  // parameters from configuration
  patElectronsInputTag_  = iConfig.getParameter<edm::InputTag>("patElectronsInputTag");
  recoElectronsInputTag_ = iConfig.getParameter<edm::InputTag>("recoElectronsInputTag");

}


PATElectronMaker::~PATElectronMaker() {}
 

// ------------ method called to produce the data  ------------
void PATElectronMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // get jet collection
  edm::Handle<vector<pat::Electron> > patElectronsHandle;
  iEvent.getByLabel(patElectronsInputTag_, patElectronsHandle);
  vector<pat::Electron> v_patElectrons = *(patElectronsHandle.product());

  edm::Handle<vector<reco::GsfElectron> > recoElectronsHandle;
  iEvent.getByLabel(recoElectronsInputTag_, recoElectronsHandle);
  vector<reco::GsfElectron> v_recoElectrons = *(recoElectronsHandle.product());

  MatchUtilities::alignRecoPatElectronCollections(v_recoElectrons, v_patElectrons);

  // create containers
  auto_ptr<vector<LorentzVector>  >   els_pat_p4               (new vector<LorentzVector>     );
  auto_ptr<vector<int>            >   els_pat_genID            (new vector<int>               );
  auto_ptr<vector<int>            >   els_pat_genMotherID      (new vector<int>               );
  auto_ptr<vector<uint32_t>       >   els_pat_flag             (new vector<uint32_t>          );
  auto_ptr<vector<float>          >   els_pat_trackIso         (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_caloIso          (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_ecalIso          (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_hcalIso          (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_robustLooseId    (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_robustTightId    (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_robustHighEnergy (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_looseId          (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_tightId          (new vector<float>             );
  auto_ptr<vector<LorentzVector>  >   els_pat_genP4            (new vector<LorentzVector>     );
  auto_ptr<vector<LorentzVector>  >   els_pat_genMotherP4      (new vector<LorentzVector>     );
  auto_ptr<vector<float>          >   els_pat_sigmaEtaEta      (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_sigmaIEtaIEta    (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_scE1x5           (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_scE2x5Max        (new vector<float>             );
  auto_ptr<vector<float>          >   els_pat_scE5x5           (new vector<float>             );
 
  // loop over top electrons and fill containers
  for ( vector<pat::Electron>::const_iterator patel_it = v_patElectrons.begin();
	patel_it != v_patElectrons.end();
	patel_it++) {
    
    els_pat_p4->push_back( patel_it->p4() );

    GenParticle gen(patel_it->genLepton() ? *patel_it->genLepton() : 
		    reco::GenParticle(0, reco::Particle::LorentzVector(0, 0, 0, 0), reco::Particle::Point(0,0,0), 0, 0, true));
    const GenParticle *gen_mom = MCUtilities::motherID(gen);
    els_pat_genID             ->push_back( gen.pdgId()                     );
    els_pat_genMotherID       ->push_back( gen_mom->pdgId()                );
    els_pat_flag              ->push_back( patel_it->status()              );
    els_pat_trackIso          ->push_back( patel_it->trackIso()            );
    els_pat_caloIso           ->push_back( patel_it->caloIso()             );
    els_pat_ecalIso           ->push_back( patel_it->ecalIso()             );
    els_pat_hcalIso           ->push_back( patel_it->hcalIso()             );
    els_pat_robustLooseId     ->push_back( patel_it->electronID("eidRobustLoose")    );
    els_pat_robustTightId     ->push_back( patel_it->electronID("eidRobustTight")    );
    els_pat_looseId           ->push_back( patel_it->electronID("eidLoose")    );
    els_pat_tightId           ->push_back( patel_it->electronID("eidTight")    );
    els_pat_robustHighEnergy  ->push_back( patel_it->electronID("eidRobustHighEnergy") );
    els_pat_genP4             ->push_back( gen.p4()                        );
    els_pat_genMotherP4       ->push_back( gen_mom->p4()                   );
    els_pat_sigmaEtaEta       ->push_back( patel_it->scSigmaEtaEta()       );
    els_pat_sigmaIEtaIEta     ->push_back( patel_it->scSigmaIEtaIEta()     );
    els_pat_scE1x5            ->push_back( patel_it->scE1x5()              );
    els_pat_scE2x5Max         ->push_back( patel_it->scE2x5Max()           );
    els_pat_scE5x5            ->push_back( patel_it->scE5x5()              );
    
    
  }

  // put containers into event
  iEvent.put(els_pat_p4,                    "elspatp4"                      );
  iEvent.put(els_pat_genID,             "elspatgenID"                );
  iEvent.put(els_pat_genMotherID,       "elspatgenMotherID"          );
  iEvent.put(els_pat_flag,              "elspatflag"                 );
  iEvent.put(els_pat_trackIso,          "elspattrackIso"             );
  iEvent.put(els_pat_caloIso,           "elspatcaloIso"              );
  iEvent.put(els_pat_ecalIso,           "elspatecalIso"              );
  iEvent.put(els_pat_hcalIso,           "elspathcalIso"              );
  iEvent.put(els_pat_robustLooseId,     "elspatrobustLooseId"        );
  iEvent.put(els_pat_robustTightId,     "elspatrobustTightId"        );
  iEvent.put(els_pat_looseId,           "elspatlooseId"              );
  iEvent.put(els_pat_tightId,           "elspattightId"              );

  iEvent.put(els_pat_robustHighEnergy,  "elspatrobustHighEnergy"     );
  iEvent.put(els_pat_genP4,             "elspatgenP4"                );
  iEvent.put(els_pat_genMotherP4,       "elspatgenMotherP4"          );
  iEvent.put(els_pat_sigmaEtaEta,       "elspatsigmaEtaEta"          );
  iEvent.put(els_pat_sigmaIEtaIEta,     "elspatsigmaIEtaIEta"        );
  iEvent.put(els_pat_scE1x5,            "elspatscE1x5"               );
  iEvent.put(els_pat_scE2x5Max,         "elspatscE2x5Max"            );
  iEvent.put(els_pat_scE5x5,            "elspatscE5x5"               );
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
PATElectronMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATElectronMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATElectronMaker);
