// -*- C++ -*-
//
// Package:    TQElectronMaker
// Class:      TQElectronMaker
// 
/**\class TQElectronMaker TQElectronMaker.cc CMS2/NtupleMaker/src/TQElectronMaker.cc

Description: copy additional TQAF jet variables in simple data structures into the EDM event tree

 Implementation:
     - take TQAF jets
     - extract and fill variables
*/
//
// Original Author:  Puneeth Kalavase
// Thu Jun 12 22:55:46 UTC 2008
// $Id: TQElectronMaker.cc,v 1.2 2008/06/14 19:05:24 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/TQElectronMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "AnalysisDataFormats/TopObjects/interface/TopElectron.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;



TQElectronMaker::TQElectronMaker(const edm::ParameterSet& iConfig) {

  // product of this EDProducer
  produces<vector<int>            >   ("elstqegammaTkNumIso"      ).setBranchAlias("els_tq_egammaTkNumIso"   );
  produces<vector<int>            >   ("elstqgenID"               ).setBranchAlias("els_tq_genID"            );
  //produces<vector<int>            >   ("elstqgenMotherID"         ).setBranchAlias("els_tq_genMotherID"      );
  produces<vector<float>          >   ("elstqtrackIso"            ).setBranchAlias("els_tq_trackIso"         );
  produces<vector<float>          >   ("elstqcaloIso"             ).setBranchAlias("els_tq_caloIso"          );
  produces<vector<float>          >   ("elstqleptonID"            ).setBranchAlias("els_tq_leptonID"         );
  produces<vector<float>          >   ("elstqelectronIDRobust"    ).setBranchAlias("els_tq_electronIDRobust" );
  produces<vector<float>          >   ("elstqegammaTkIso"         ).setBranchAlias("els_tq_egammaTkIso"      );
  produces<vector<float>          >   ("elstqegammaEcalIso"       ).setBranchAlias("els_tq_egammaEcalIso"    );
  produces<vector<float>          >   ("elstqegammaHcalIso"       ).setBranchAlias("els_tq_egammaHcalIso"    );
  produces<vector<float>          >   ("elstqLRComb"              ).setBranchAlias("els_tq_LRComb"           );
  produces<vector<LorentzVector>  >   ("elstqgenP4"               ).setBranchAlias("els_tq_genP4"            );
  //produces<vector<LorentzVector>  >   ("elstqgenMtherP4"          ).setBranchAlias("els_tq_genMotherP4"      );
  

  // parameters from configuration
  tqElectronsInputTag = iConfig.getParameter<edm::InputTag>("tqElectronsInputTag");

}


TQElectronMaker::~TQElectronMaker() {}
 

// ------------ method called to produce the data  ------------
void TQElectronMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // get jet collection
  edm::Handle<std::vector<TopElectron> > tqElectronHandle;
  iEvent.getByLabel(tqElectronsInputTag, tqElectronHandle);

  // create containers
  auto_ptr<vector<int>            >   els_tq_egammaTkNumIso(   new vector<int>               );
  auto_ptr<vector<int>            >   els_tq_genID(            new vector<int>               );
  //auto_ptr<vector<int>            >   els_tq_genMotherID(      new vector<int>               );
   auto_ptr<vector<float>          >   els_tq_trackIso(         new vector<float>             );
   auto_ptr<vector<float>          >   els_tq_caloIso(          new vector<float>             );
   auto_ptr<vector<float>          >   els_tq_leptonID(         new vector<float>             );
   auto_ptr<vector<float>          >   els_tq_electronIDRobust( new vector<float>             );
   auto_ptr<vector<float>          >   els_tq_egammaTkIso(      new vector<float>             );
   auto_ptr<vector<float>          >   els_tq_egammaEcalIso(    new vector<float>             );
   auto_ptr<vector<float>          >   els_tq_egammaHcalIso(    new vector<float>             );
   auto_ptr<vector<float>          >   els_tq_LRComb(           new vector<float>             ); 
   auto_ptr<vector<LorentzVector>  >   els_tq_genP4(            new vector<LorentzVector>     );
  //auto_ptr<vector<LorentzVector>  >   els_tq_genMotherP4(      new vector<LorentzVector>     );
  
  // loop over top electrons and fill containers
  for ( std::vector<TopElectron>::const_iterator tqel_it = tqElectronHandle->begin();
	tqel_it != tqElectronHandle->end();
	tqel_it++) {
    
    GenParticle gen(tqel_it->getGenLepton());
    
    els_tq_egammaTkNumIso    ->push_back( tqel_it->getEgammaTkNumIso()   );
    els_tq_genID             ->push_back( gen.pdgId()                    );
    //els_tq_genMotherID       ->push_back( gen.
    els_tq_trackIso          ->push_back( tqel_it->getTrackIso()         );
    els_tq_caloIso           ->push_back( tqel_it->getCaloIso()          );
    els_tq_leptonID          ->push_back( tqel_it->getLeptonID()         );
    els_tq_electronIDRobust  ->push_back( tqel_it->getElectronIDRobust() );
    els_tq_egammaTkIso       ->push_back( tqel_it->getEgammaTkIso()      );
    els_tq_egammaEcalIso     ->push_back( tqel_it->getEgammaEcalIso()    );
     els_tq_egammaHcalIso     ->push_back( tqel_it->getEgammaHcalIso()    );
     els_tq_LRComb            ->push_back( tqel_it->getLRComb()           );
     els_tq_genP4             ->push_back( gen.p4()                       );
    //els_tq_genMotherP4       ->push_back( 
    
  }

  // put containers into event
  iEvent.put(els_tq_egammaTkNumIso,    "elstqegammaTkNumIso"       );
  iEvent.put(els_tq_genID,             "elstqgenID"                );
  //iEvent.put(els_tq_genMotherID,       "elstqgenMotherID"          );
  iEvent.put(els_tq_trackIso,          "elstqtrackIso"             );
  iEvent.put(els_tq_caloIso,           "elstqcaloIso"              );
  iEvent.put(els_tq_leptonID,          "elstqleptonID"             );
  iEvent.put(els_tq_electronIDRobust,  "elstqelectronIDRobust"     );
  iEvent.put(els_tq_egammaTkIso,       "elstqegammaTkIso"          );
  iEvent.put(els_tq_egammaEcalIso,     "elstqegammaEcalIso"        );
  iEvent.put(els_tq_egammaHcalIso,     "elstqegammaHcalIso"        );
  iEvent.put(els_tq_LRComb,            "elstqLRComb"               );
  iEvent.put(els_tq_genP4,             "elstqgenP4"                );
  //iEvent.put(els_tq_genMotherP4,       "elstqgenMotherP4"          );
}

// ------------ method called once each job just before starting event loop  ------------
void 
TQElectronMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TQElectronMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TQElectronMaker);
