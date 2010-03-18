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
// $Id: PATElectronMaker.cc,v 1.9 2010/03/18 02:13:15 kalavase Exp $
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

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;



PATElectronMaker::PATElectronMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // product of this EDProducer
  produces<vector<LorentzVector>  >   (branchprefix+"p4").setBranchAlias(aliasprefix_+"_p4");
  produces<vector<int>            >   (branchprefix+"genID"               ).setBranchAlias(aliasprefix_+"_genID"            );
  produces<vector<int>            >   (branchprefix+"genMotherID"         ).setBranchAlias(aliasprefix_+"_genMotherID"      );
  produces<vector<uint32_t>       >   (branchprefix+"flag"                ).setBranchAlias(aliasprefix_+"_flag"             );
  produces<vector<float>          >   (branchprefix+"trackIso"            ).setBranchAlias(aliasprefix_+"_trackIso"         );
  produces<vector<float>          >   (branchprefix+"caloIso"             ).setBranchAlias(aliasprefix_+"_caloIso"          );
  produces<vector<float>          >   (branchprefix+"ecalIso"             ).setBranchAlias(aliasprefix_+"_ecalIso"          );
  produces<vector<float>          >   (branchprefix+"hcalIso"             ).setBranchAlias(aliasprefix_+"_hcalIso"          );
  produces<vector<float>          >   (branchprefix+"robustLooseId"       ).setBranchAlias(aliasprefix_+"_robustLooseId"    );
  produces<vector<float>          >   (branchprefix+"robustTightId"       ).setBranchAlias(aliasprefix_+"_robustTightId"    );
  produces<vector<float>          >   (branchprefix+"looseId"             ).setBranchAlias(aliasprefix_+"_looseId"          );
  produces<vector<float>          >   (branchprefix+"tightId"             ).setBranchAlias(aliasprefix_+"_tightId"          );
  produces<vector<float>          >   (branchprefix+"robustHighEnergy"    ).setBranchAlias(aliasprefix_+"_robustHighEnergy" );
  produces<vector<LorentzVector>  >   (branchprefix+"genP4"               ).setBranchAlias(aliasprefix_+"_genP4"            );
  produces<vector<LorentzVector>  >   (branchprefix+"genMotherP4"         ).setBranchAlias(aliasprefix_+"_genMotherP4"      );
  produces<vector<float>          >   (branchprefix+"sigmaEtaEta"         ).setBranchAlias(aliasprefix_+"_sigmaEtaEta"      );
  produces<vector<float>          >   (branchprefix+"sigmaIEtaIEta"       ).setBranchAlias(aliasprefix_+"_sigmaIEtaIEta"    );
  produces<vector<float>          >   (branchprefix+"scE1x5"              ).setBranchAlias(aliasprefix_+"_scE1x5"           );
  produces<vector<float>          >   (branchprefix+"scE2x5Max"           ).setBranchAlias(aliasprefix_+"_scE2x5Max"        );
  produces<vector<float>          >   (branchprefix+"scE5x5"              ).setBranchAlias(aliasprefix_+"_scE5x5"           );

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
    
    els_pat_p4->push_back( LorentzVector( patel_it->p4() ) );

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
    els_pat_genP4             ->push_back( LorentzVector( gen.p4() )       );
    els_pat_genMotherP4       ->push_back( LorentzVector( gen_mom->p4() )  );
    els_pat_sigmaEtaEta       ->push_back( patel_it->scSigmaEtaEta()       );
    els_pat_sigmaIEtaIEta     ->push_back( patel_it->scSigmaIEtaIEta()     );
    els_pat_scE1x5            ->push_back( patel_it->scE1x5()              );
    els_pat_scE2x5Max         ->push_back( patel_it->scE2x5Max()           );
    els_pat_scE5x5            ->push_back( patel_it->scE5x5()              );
    
    
  }

  // put containers into event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(els_pat_p4,                    branchprefix+"p4"                      );
  iEvent.put(els_pat_genID,             branchprefix+"genID"                );
  iEvent.put(els_pat_genMotherID,       branchprefix+"genMotherID"          );
  iEvent.put(els_pat_flag,              branchprefix+"flag"                 );
  iEvent.put(els_pat_trackIso,          branchprefix+"trackIso"             );
  iEvent.put(els_pat_caloIso,           branchprefix+"caloIso"              );
  iEvent.put(els_pat_ecalIso,           branchprefix+"ecalIso"              );
  iEvent.put(els_pat_hcalIso,           branchprefix+"hcalIso"              );
  iEvent.put(els_pat_robustLooseId,     branchprefix+"robustLooseId"        );
  iEvent.put(els_pat_robustTightId,     branchprefix+"robustTightId"        );
  iEvent.put(els_pat_looseId,           branchprefix+"looseId"              );
  iEvent.put(els_pat_tightId,           branchprefix+"tightId"              );

  iEvent.put(els_pat_robustHighEnergy,  branchprefix+"robustHighEnergy"     );
  iEvent.put(els_pat_genP4,             branchprefix+"genP4"                );
  iEvent.put(els_pat_genMotherP4,       branchprefix+"genMotherP4"          );
  iEvent.put(els_pat_sigmaEtaEta,       branchprefix+"sigmaEtaEta"          );
  iEvent.put(els_pat_sigmaIEtaIEta,     branchprefix+"sigmaIEtaIEta"        );
  iEvent.put(els_pat_scE1x5,            branchprefix+"scE1x5"               );
  iEvent.put(els_pat_scE2x5Max,         branchprefix+"scE2x5Max"            );
  iEvent.put(els_pat_scE5x5,            branchprefix+"scE5x5"               );
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
PATElectronMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATElectronMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATElectronMaker);
