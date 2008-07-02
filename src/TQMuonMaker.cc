// -*- C++ -*-
//
// Package:    TQMuonMaker
// Class:      TQMuonMaker
// 
/**\class TQMuonMaker TQMuonMaker.cc CMS2/NtupleMaker/src/TQMuonMaker.cc

Description: copy additional TQAF muon variables in simple data structures into the EDM event tree

 Implementation:
     - take TQAF muons
     - extract and fill variables
*/
//
// Original Author:  Frank Golf
// Thu Jun 25 16:39:55 UTC 2008
// $Id: TQMuonMaker.cc,v 1.1 2008/07/02 02:26:12 jmuelmen Exp $
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
#include "CMS2/NtupleMaker/interface/TQMuonMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "AnalysisDataFormats/TopObjects/interface/TopMuon.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;

TQMuonMaker::TQMuonMaker(const edm::ParameterSet& iConfig) {

  // product of this EDProducer
  produces<vector<float>          >   ("mustqtrackIso"    ).setBranchAlias("mus_tq_trackIso"    );
  produces<vector<float>          >   ("mustqcaloIso"     ).setBranchAlias("mus_tq_caloIso"     );
  produces<vector<float>          >   ("mustqleptonID"    ).setBranchAlias("mus_tq_leptonID"    );
  produces<vector<float>          >   ("mustqlrComb"      ).setBranchAlias("mus_tq_lrComb"      );
  produces<vector<int>            >   ("mustqgenID"       ).setBranchAlias("mus_tq_genID"       );
  produces<vector<int>            >   ("mustqgenMotherID" ).setBranchAlias("mus_tq_genMotherID" );
  produces<vector<LorentzVector>  >   ("mustqgenP4"       ).setBranchAlias("mus_tq_genP4"       );
  produces<vector<LorentzVector>  >   ("mustqgenMotherP4" ).setBranchAlias("mus_tq_genMotherP4" );

  // parameters from configuration
  tqMuonsInputTag = iConfig.getParameter<edm::InputTag>("tqMuonsInputTag");
}

TQMuonMaker::~TQMuonMaker() {}

// ------------ method called to produce the data  ------------
void TQMuonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // get jet collection
  edm::Handle<std::vector<TopMuon> > tqMuonHandle;
  iEvent.getByLabel(tqMuonsInputTag, tqMuonHandle);

  // create containers
  auto_ptr<vector<float>          >   mus_tq_trackIso    ( new vector<float>             );
  auto_ptr<vector<float>          >   mus_tq_caloIso     ( new vector<float>             );
  auto_ptr<vector<float>          >   mus_tq_leptonID    ( new vector<float>             );
  auto_ptr<vector<float>          >   mus_tq_lrComb      ( new vector<float>             );
  auto_ptr<vector<int>            >   mus_tq_genID       ( new vector<int>               );
  auto_ptr<vector<int>            >   mus_tq_genMotherID ( new vector<int>               ); 
  auto_ptr<vector<LorentzVector>  >   mus_tq_genP4       ( new vector<LorentzVector>     );
  auto_ptr<vector<LorentzVector>  >   mus_tq_genMotherP4 ( new vector<LorentzVector>     );
  
  // loop over top muons and fill containers
  for ( std::vector<TopMuon>::const_iterator tqmu_it = tqMuonHandle->begin(),
	tqmu_end = tqMuonHandle->end();
	tqmu_it != tqmu_end; tqmu_it++) {
    
    GenParticle gen( tqmu_it->getGenLepton() );
    const GenParticle *gen_mom = MCUtilities::motherID(gen);
    mus_tq_trackIso    ->push_back( tqmu_it->getTrackIso() );
    mus_tq_caloIso     ->push_back( tqmu_it->getCaloIso() );
    mus_tq_leptonID    ->push_back( tqmu_it->getLeptonID() );
    mus_tq_lrComb      ->push_back( tqmu_it->getLRComb() );
    mus_tq_genID       ->push_back( gen.pdgId() );
    mus_tq_genMotherID ->push_back( gen_mom->pdgId() );
    mus_tq_genP4       ->push_back( gen.p4() );
    mus_tq_genMotherP4 ->push_back( gen_mom->p4() );
  }

  // put containers into event
  iEvent.put(mus_tq_trackIso,     "mustqtrackIso"        );
  iEvent.put(mus_tq_caloIso,      "mustqcaloIso"         );
  iEvent.put(mus_tq_leptonID,     "mustqleptonID"        );
  iEvent.put(mus_tq_lrComb,       "mustqlrComb"          );
  iEvent.put(mus_tq_genID,        "mustqgenID"           );
  iEvent.put(mus_tq_genMotherID,  "mustqgenMotherID"     );
  iEvent.put(mus_tq_genP4,        "mustqgenP4"           );
  iEvent.put(mus_tq_genMotherP4,  "mustqgenMotherP4"     );
}

// ------------ method called once each job just before starting event loop  ------------
void 
TQMuonMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TQMuonMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TQMuonMaker);
