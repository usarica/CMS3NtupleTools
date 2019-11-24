// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      HypDilepMaker
//
// Original Author:  Puneeth Kalavase
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
/*
Description: create trilepton hypothesis branches

Implementation:
- combine muons and electrons after preselection
- correct met

Hypothesis (lt):
mm:0
me:1 (only happens if m is > tight cut and e < tight cut. Avoids double counting) 
em:2 
ee:3

*/

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS3/NtupleMaker/interface/plugins/HypDilepMaker.h"
#include "CMS3/NtupleMaker/interface/plugins/MatchUtilities.h"

#include "Math/VectorUtil.h"
#include "TMath.h"


using namespace reco;
using namespace edm;
using namespace std;


HypDilepMaker::HypDilepMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
  
  musChargeToken = consumes<std::vector<int> >(iConfig.getParameter<InputTag>("musChargeInputTag"    ));
  musTypeToken   = consumes<std::vector<int> >(iConfig.getParameter<InputTag>("musTypeInputTag"    ));
  musp4Token     = consumes<std::vector<LorentzVector> >(iConfig.getParameter<InputTag>("musp4InputTag"    ));

  elsChargeToken = consumes<std::vector<int> >(iConfig.getParameter<InputTag>("elsChargeInputTag"    ));
  elsTypeToken   = consumes<std::vector<int> >(iConfig.getParameter<InputTag>("elsTypeInputTag"    ));
  elsp4Token     = consumes<std::vector<LorentzVector> >(iConfig.getParameter<InputTag>("elsp4InputTag"    ));

  tightptcut               = iConfig.getParameter<double>  ("TightLepton_PtCut");
  looseptcut               = iConfig.getParameter<double>  ("LooseLepton_PtCut");

  produces<vector<int> >           (branchprefix+"type"    ).setBranchAlias(aliasprefix_+"_type"     );
  produces<vector<LorentzVector> > (branchprefix+"p4"      ).setBranchAlias(aliasprefix_+"_p4"       );
  
  produces<vector<int> >           (branchprefix+"ltcharge").setBranchAlias(aliasprefix_+"_lt_charge");
  produces<vector<int> >           (branchprefix+"ltindex" ).setBranchAlias(aliasprefix_+"_lt_index" );
  produces<vector<int> >           (branchprefix+"ltid"    ).setBranchAlias(aliasprefix_+"_lt_id"    );
  produces<vector<LorentzVector > >(branchprefix+"ltp4"    ).setBranchAlias(aliasprefix_+"_lt_p4"    );
  
  produces<vector<int> >           (branchprefix+"llcharge").setBranchAlias(aliasprefix_+"_ll_charge");
  produces<vector<int> >           (branchprefix+"llindex" ).setBranchAlias(aliasprefix_+"_ll_index" );
  produces<vector<int> >           (branchprefix+"llid"    ).setBranchAlias(aliasprefix_+"_ll_id"    );
  produces<vector<LorentzVector > >(branchprefix+"llp4"    ).setBranchAlias(aliasprefix_+"_ll_p4"    );
  
}

HypDilepMaker::~HypDilepMaker() {}

// ------------ method called to produce the data  ------------
void HypDilepMaker::produce(Event& iEvent, const edm::EventSetup& iSetup) {

  // output collections
  unique_ptr<vector<int> >           hyp_type                    (new vector<int>             );
  unique_ptr<vector<LorentzVector> > hyp_p4                      (new vector<LorentzVector>   );
  unique_ptr<vector<int> >           hyp_lt_charge               (new vector<int>             );
  unique_ptr<vector<int> >           hyp_lt_index                (new vector<int>             );
  unique_ptr<vector<int> >           hyp_lt_id                   (new vector<int>             );
  unique_ptr<vector<LorentzVector> > hyp_lt_p4                   (new vector<LorentzVector>   );
  unique_ptr<vector<int> >           hyp_ll_charge               (new vector<int>             );
  unique_ptr<vector<int> >           hyp_ll_index                (new vector<int>             );
  unique_ptr<vector<int> >           hyp_ll_id                   (new vector<int>             );
  unique_ptr<vector<LorentzVector> > hyp_ll_p4                   (new vector<LorentzVector>   );
  
  // muon charge
  // edm::InputTag mus_charge_tag(muonsInputTag.label(),"muscharge");
  edm::Handle<std::vector<int> > mus_charge_h;
  iEvent.getByToken(musChargeToken, mus_charge_h);
  const vector<int> *mus_charge = mus_charge_h.product();

  //muon p4
  // InputTag mus_p4_tag(muonsInputTag.label(),"musp4");
  Handle<vector<LorentzVector> > mus_p4_h;
  iEvent.getByToken(musp4Token, mus_p4_h);
  const vector<LorentzVector> *mus_p4 = mus_p4_h.product();

  //muon type
  // InputTag mus_type_tag(muonsInputTag.label(), "mustype");
  Handle<vector<int> > mus_type_h;
  iEvent.getByToken(musTypeToken, mus_type_h);
  const vector<int> *mus_type = mus_type_h.product();

  //-----------------------------------------------------------
  // electron variables
  //-----------------------------------------------------------
  // InputTag els_charge_tag(electronsInputTag.label(),"elscharge");
  Handle<vector<int> > els_charge_h;
  iEvent.getByToken(elsChargeToken, els_charge_h);
  const vector<int> *els_charge = els_charge_h.product();

  // electron p4
  // InputTag els_p4_tag(electronsInputTag.label(),"elsp4");
  Handle<vector<LorentzVector> > els_p4_h;
  iEvent.getByToken(elsp4Token, els_p4_h);
  const vector<LorentzVector> *els_p4 = els_p4_h.product();

  unsigned int nmus = mus_p4->size();
  unsigned int nels = els_p4->size();

  //------------------------------------------------------------
  // loop over the muons
  //------------------------------------------------------------
  //get the candidates and make hypotheses 
  for(unsigned int mus_index_1 = 0; mus_index_1 < nmus; mus_index_1++) {//first muon loop
    for(unsigned int mus_index_2 = 0; mus_index_2 < nmus; mus_index_2++) {//second muon loop

      if(mus_index_1 == mus_index_2) continue;
      if(mus_index_2 < mus_index_1)  continue;  //avoid double counting

      //don't look at standalone muons
      if(mus_type->at(mus_index_1) == 8) continue;
      if(mus_type->at(mus_index_2) == 8) continue;
      
      float mu_pt1 = mus_p4->at(mus_index_1).Pt();
      float mu_pt2 = mus_p4->at(mus_index_2).Pt();
      
      //if either fail the loose cut, go to the next muon
      if(mu_pt1 < looseptcut || mu_pt2 < looseptcut) continue;
      
      //if neither one passes the tight cut, go to the next muon
      if(mu_pt1 < tightptcut && mu_pt2 < tightptcut) continue;

      int tight_index = mus_index_1;
      int loose_index = mus_index_2;

      /*
	figure out which one should be tight and which should
	be loose in case one passes the tight cut and the other 
	does not
      */
      if(mu_pt1 < tightptcut && mu_pt2 > tightptcut) {
	tight_index = mus_index_2;
	loose_index = mus_index_1;
      }
      if(mu_pt2 < tightptcut && mu_pt1 > tightptcut) {
	tight_index = mus_index_1;
	loose_index = mus_index_2;
      }
            
      hyp_type          ->push_back(0                                     );
      hyp_p4            ->push_back(mus_p4->at(tight_index)+mus_p4->at(loose_index)               );
      hyp_lt_charge       ->push_back(mus_charge       ->at(tight_index)  );
      hyp_lt_index        ->push_back(tight_index                         );
      hyp_lt_id           ->push_back(-13*(mus_charge   ->at(tight_index)));
      hyp_lt_p4           ->push_back(mus_p4           ->at(tight_index)  );
      hyp_ll_charge       ->push_back(mus_charge       ->at(loose_index)  );
      hyp_ll_index        ->push_back(loose_index                         );
      hyp_ll_id           ->push_back(-13*(mus_charge   ->at(loose_index)));
      hyp_ll_p4           ->push_back(mus_p4           ->at(loose_index)  );
    }
  }  

  //------------------------------------------------------------
  // loop over the elecrons
  //------------------------------------------------------------
  //get the candidates and make hypotheses 
  for(unsigned int els_index_1 = 0; els_index_1 < nels; els_index_1++) {
    for(unsigned int els_index_2 = 0; els_index_2 < nels; els_index_2++) {
      
      if(els_index_1 == els_index_2) continue;
      if(els_index_2 < els_index_1)  continue;  //avoid double counting
      
      float el_pt1 = els_p4->at(els_index_1).Pt();
      float el_pt2 = els_p4->at(els_index_2).Pt();
      
      //if either fail the loose cut, go to the next muon
      if(el_pt1 < looseptcut || el_pt2 < looseptcut) continue;
      
      //if neither one passes the tight cut, continue
      if(el_pt1 < tightptcut && el_pt2 < tightptcut) continue;
      
      int tight_index = els_index_1;
      int loose_index = els_index_2;
      
      /*
	figure out which one should be tight and which should
	be loose in case one passes the tight cut and the other 
	does not
      */
      if(el_pt1 < tightptcut && el_pt2 > tightptcut) {
	tight_index = els_index_2;
	loose_index = els_index_1;
      }
      if(el_pt2 < tightptcut && el_pt1 > tightptcut) {
	tight_index = els_index_1;
	loose_index = els_index_2;
      }

      //hyp_Ht->push_back(temp_Ht);
    
      hyp_type          ->push_back(3);
      hyp_p4            ->push_back(els_p4->at(tight_index)+els_p4->at(loose_index)               );
      hyp_lt_charge       ->push_back(els_charge       ->at(tight_index)  );
      hyp_lt_index        ->push_back(tight_index                         );
      hyp_lt_id           ->push_back(-11*(els_charge   ->at(tight_index)));
      hyp_lt_p4           ->push_back(els_p4           ->at(tight_index)  );
      hyp_ll_charge       ->push_back(els_charge       ->at(loose_index)  );
      hyp_ll_index        ->push_back(loose_index                         );
      hyp_ll_id           ->push_back(-11*(els_charge   ->at(loose_index)));
      hyp_ll_p4           ->push_back(els_p4           ->at(loose_index)  );
    }
  }  
  
  /*------------------------------------------------------------
    The EMu, MuE cases
    To avoid double counting, only make MuE if Mu is tight and E is loose
  */

  for(unsigned int els_index = 0; els_index < nels; els_index++) {
    for(unsigned int mus_index = 0; mus_index < nmus; mus_index++) {

      if(mus_type->at(mus_index) == 8) continue;

      float el_pt = els_p4->at(els_index).Pt();
      float mu_pt = mus_p4->at(mus_index).Pt();

      //if either fail the loose cut, go to the next muon
      if(el_pt < looseptcut || mu_pt < looseptcut) continue;

      //if both fail the tight cut, continue
      if(el_pt < tightptcut && mu_pt < tightptcut) continue;
                  
      hyp_p4->push_back(mus_p4->at(mus_index)+els_p4->at(els_index));
	
      if(el_pt < tightptcut && mu_pt > tightptcut) {
	hyp_type            ->push_back(1);
	  
	hyp_lt_charge       ->push_back(mus_charge       ->at(mus_index)  );
	hyp_lt_index        ->push_back(mus_index                         );
	hyp_lt_id           ->push_back(-13*(mus_charge   ->at(mus_index)));
	hyp_lt_p4           ->push_back(mus_p4           ->at(mus_index)  );
	hyp_ll_charge       ->push_back(els_charge       ->at(els_index)  );
	hyp_ll_index        ->push_back(els_index                         );
	hyp_ll_id           ->push_back(-11*(els_charge   ->at(els_index)));
	hyp_ll_p4           ->push_back(els_p4           ->at(els_index)  );
      }
    else {
	hyp_type            ->push_back(2);
	hyp_lt_charge       ->push_back(els_charge       ->at(els_index)  );
	hyp_lt_index        ->push_back(els_index                         );
	hyp_lt_id           ->push_back(-11*(els_charge   ->at(els_index)));
	hyp_lt_p4           ->push_back(els_p4           ->at(els_index)  );
	
	hyp_ll_charge       ->push_back(mus_charge       ->at(mus_index)  );
	hyp_ll_index        ->push_back(mus_index                         );
	hyp_ll_id           ->push_back(-13*(mus_charge   ->at(mus_index)));
	hyp_ll_p4           ->push_back(mus_p4           ->at(mus_index)  );
	
      }
    }
  }

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(std::move(hyp_type                     ),branchprefix+"type"                     );
  iEvent.put(std::move(hyp_p4                       ),branchprefix+"p4"                       );
  iEvent.put(std::move(hyp_lt_charge                ),branchprefix+"ltcharge"                 );
  iEvent.put(std::move(hyp_lt_index                 ),branchprefix+"ltindex"                  );
  iEvent.put(std::move(hyp_lt_id                    ),branchprefix+"ltid"                     );
  iEvent.put(std::move(hyp_lt_p4                    ),branchprefix+"ltp4"                     );
  iEvent.put(std::move(hyp_ll_charge                ),branchprefix+"llcharge"                 );
  iEvent.put(std::move(hyp_ll_index                 ),branchprefix+"llindex"                  );
  iEvent.put(std::move(hyp_ll_id                    ),branchprefix+"llid"                     );
  iEvent.put(std::move(hyp_ll_p4                    ),branchprefix+"llp4"                     );
}

// ------------ method called once each job just before starting event loop  ------------
void HypDilepMaker::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void HypDilepMaker::endJob(){
}

//define this as a plug-in
DEFINE_FWK_MODULE(HypDilepMaker);
