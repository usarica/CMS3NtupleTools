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
- correct jets and store index vectors
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
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS3/NtupleMaker/interface/HypDilepMaker.h"
#include "CMS3/NtupleMaker/interface/MatchUtilities.h"

#include "CMS3/NtupleMaker/interface/MT2Utility.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "Math/VectorUtil.h"
#include "TMath.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;


HypDilepMaker::HypDilepMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
  
  muonsInputTag            = iConfig.getParameter<InputTag>("muonsInputTag"                                        );
  electronsInputTag        = iConfig.getParameter<InputTag>("electronsInputTag"                                    );
  metInputTag              = iConfig.getParameter<InputTag>("metInputTag"                                          );
  jetsInputTag             = iConfig.getParameter<InputTag>("jetsInputTag"                                         );
  hypJetMaxEtaCut          = iConfig.getParameter<double>  ("hypJetMaxEtaCut"                                      );
  hypJetMinPtCut           = iConfig.getParameter<double>  ("hypJetMinPtCut"                                       );
  tightptcut               = iConfig.getParameter<double>  ("TightLepton_PtCut"                                    );
  looseptcut               = iConfig.getParameter<double>  ("LooseLepton_PtCut"                                    );

  produces<vector<int> >           (branchprefix+"type"                    ).setBranchAlias(aliasprefix_+"_type"                       );
  produces<vector<LorentzVector> > (branchprefix+"p4"                      ).setBranchAlias(aliasprefix_+"_p4"                         );
  
  produces<vector<int> >           (branchprefix+"ltcharge"                ).setBranchAlias(aliasprefix_+"_lt_charge"                  );
  produces<vector<int> >           (branchprefix+"ltindex"                 ).setBranchAlias(aliasprefix_+"_lt_index"                   );
  produces<vector<int> >           (branchprefix+"ltid"                    ).setBranchAlias(aliasprefix_+"_lt_id"                      );
  produces<vector<float> >         (branchprefix+"ltd0"                    ).setBranchAlias(aliasprefix_+"_lt_d0"                      );
  produces<vector<float> >         (branchprefix+"ltz0"                    ).setBranchAlias(aliasprefix_+"_lt_z0"                      );
  produces<vector<float> >         (branchprefix+"ltd0corr"                ).setBranchAlias(aliasprefix_+"_lt_d0corr"                  );
  produces<vector<float> >         (branchprefix+"ltz0corr"                ).setBranchAlias(aliasprefix_+"_lt_z0corr"                  );
  produces<vector<float> >         (branchprefix+"ltchi2"                  ).setBranchAlias(aliasprefix_+"_lt_chi2"                    );
  produces<vector<float> >         (branchprefix+"ltndof"                  ).setBranchAlias(aliasprefix_+"_lt_ndof"                    );
  produces<vector<float> >         (branchprefix+"ltd0Err"                 ).setBranchAlias(aliasprefix_+"_lt_d0Err"                   );
  produces<vector<float> >         (branchprefix+"ltz0Err"                 ).setBranchAlias(aliasprefix_+"_lt_z0Err"                   );
  produces<vector<LorentzVector > >(branchprefix+"ltp4"                    ).setBranchAlias(aliasprefix_+"_lt_p4"                      );
  produces<vector<LorentzVector > >(branchprefix+"lttrkp4"                 ).setBranchAlias(aliasprefix_+"_lt_trk_p4"                  );
  
  produces<vector<int> >           (branchprefix+"llcharge"                ).setBranchAlias(aliasprefix_+"_ll_charge"                  );
  produces<vector<int> >           (branchprefix+"llindex"                 ).setBranchAlias(aliasprefix_+"_ll_index"                   );
  produces<vector<int> >           (branchprefix+"llid"                    ).setBranchAlias(aliasprefix_+"_ll_id"                      );
  produces<vector<float> >         (branchprefix+"lld0"                    ).setBranchAlias(aliasprefix_+"_ll_d0"                      );
  produces<vector<float> >         (branchprefix+"llz0"                    ).setBranchAlias(aliasprefix_+"_ll_z0"                      );
  produces<vector<float> >         (branchprefix+"lld0corr"                ).setBranchAlias(aliasprefix_+"_ll_d0corr"                  );
  produces<vector<float> >         (branchprefix+"llz0corr"                ).setBranchAlias(aliasprefix_+"_ll_z0corr"                  );
  produces<vector<float> >         (branchprefix+"llchi2"                  ).setBranchAlias(aliasprefix_+"_ll_chi2"                    );
  produces<vector<float> >         (branchprefix+"llndof"                  ).setBranchAlias(aliasprefix_+"_ll_ndof"                    );
  produces<vector<float> >         (branchprefix+"lld0Err"                 ).setBranchAlias(aliasprefix_+"_ll_d0Err"                   );
  produces<vector<float> >         (branchprefix+"llz0Err"                 ).setBranchAlias(aliasprefix_+"_ll_z0Err"                   );
  produces<vector<LorentzVector > >(branchprefix+"llp4"                    ).setBranchAlias(aliasprefix_+"_ll_p4"                      );
  produces<vector<LorentzVector > >(branchprefix+"lltrkp4"                 ).setBranchAlias(aliasprefix_+"_ll_trk_p4"                  );
  
}


HypDilepMaker::~HypDilepMaker() {}


//
// member functions
//

// ------------ method called to produce the data  ------------
void HypDilepMaker::produce(Event& iEvent, const edm::EventSetup& iSetup) {

  // output collections
  auto_ptr<vector<int> >           hyp_type                    (new vector<int>             );
  auto_ptr<vector<LorentzVector> > hyp_p4                      (new vector<LorentzVector>   );
  auto_ptr<vector<int> >           hyp_lt_charge               (new vector<int>             );
  auto_ptr<vector<int> >           hyp_lt_index                (new vector<int>             );
  auto_ptr<vector<int> >           hyp_lt_id                   (new vector<int>             );
  auto_ptr<vector<float> >         hyp_lt_d0                   (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_z0                   (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_d0corr               (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_z0corr               (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_chi2                 (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_ndof                 (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_d0Err                (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_z0Err                (new vector<float>           );
  auto_ptr<vector<LorentzVector> > hyp_lt_p4                   (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > hyp_lt_trk_p4               (new vector<LorentzVector>   );
  auto_ptr<vector<int> >           hyp_ll_charge               (new vector<int>             );
  auto_ptr<vector<int> >           hyp_ll_index                (new vector<int>             );
  auto_ptr<vector<int> >           hyp_ll_id                   (new vector<int>             );
  auto_ptr<vector<float> >         hyp_ll_d0                   (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_z0                   (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_d0corr               (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_z0corr               (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_chi2                 (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_ndof                 (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_d0Err                (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_z0Err                (new vector<float>           );
  auto_ptr<vector<LorentzVector> > hyp_ll_p4                   (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > hyp_ll_trk_p4               (new vector<LorentzVector>   );
  
  // muon charge
  edm::InputTag mus_charge_tag(muonsInputTag.label(),"muscharge");
  edm::Handle<std::vector<int> > mus_charge_h;
  iEvent.getByLabel(mus_charge_tag, mus_charge_h);
  const vector<int> *mus_charge = mus_charge_h.product();

  //muon p4
  InputTag mus_p4_tag(muonsInputTag.label(),"musp4");
  Handle<vector<LorentzVector> > mus_p4_h;
  iEvent.getByLabel(mus_p4_tag, mus_p4_h);
  const vector<LorentzVector> *mus_p4 = mus_p4_h.product();

  //muond0
  InputTag mus_d0_tag(muonsInputTag.label(),"musd0");
  Handle<vector<float> > mus_d0_h;
  iEvent.getByLabel(mus_d0_tag, mus_d0_h);
  const vector<float> *mus_d0 = mus_d0_h.product();

  //muon d0, corrected for the beamspot
  InputTag mus_d0corr_tag(muonsInputTag.label(),"musd0corr");
  Handle<vector<float> > mus_d0corr_h;
  iEvent.getByLabel(mus_d0corr_tag, mus_d0corr_h);
  const vector<float> *mus_d0corr = mus_d0corr_h.product();

  //muon z0
  InputTag mus_z0_tag(muonsInputTag.label(),"musz0");
  Handle<vector<float> > mus_z0_h;
  iEvent.getByLabel(mus_z0_tag, mus_z0_h);
  const vector<float> *mus_z0 = mus_z0_h.product();

  //muon z0, corrected for the beamspot
  InputTag mus_z0corr_tag(muonsInputTag.label(),"musz0corr");
  Handle<vector<float> > mus_z0corr_h;
  iEvent.getByLabel(mus_z0corr_tag, mus_z0corr_h);
  const vector<float> *mus_z0corr = mus_z0corr_h.product();

  
  //chi2
  InputTag mus_chi2_tag(muonsInputTag.label(),"muschi2");
  Handle<vector<float> > mus_chi2_h;
  iEvent.getByLabel(mus_chi2_tag, mus_chi2_h);
  const vector<float> *mus_chi2 = mus_chi2_h.product();
  
  //ndof
  InputTag mus_ndof_tag(muonsInputTag.label(),"musndof");
  Handle<vector<float> > mus_ndof_h;
  iEvent.getByLabel(mus_ndof_tag, mus_ndof_h);
  const vector<float> *mus_ndof = mus_ndof_h.product();
  
  //d0err
  InputTag mus_d0err_tag(muonsInputTag.label(),"musd0Err");
  Handle<vector<float> > mus_d0err_h;
  iEvent.getByLabel(mus_d0err_tag, mus_d0err_h);
  const vector<float> *mus_d0Err = mus_d0err_h.product();

  //z0err
  InputTag mus_z0err_tag(muonsInputTag.label(),"musz0Err");
  Handle<vector<float> > mus_z0err_h;
  iEvent.getByLabel(mus_z0err_tag, mus_z0err_h);
  const vector<float> *mus_z0Err = mus_z0err_h.product();

  //muon track P4
  InputTag mus_trk_p4_tag(muonsInputTag.label(),"mustrkp4");
  Handle<vector<LorentzVector> > mus_trk_p4_h;
  iEvent.getByLabel(mus_trk_p4_tag, mus_trk_p4_h);
  const vector<LorentzVector> *mus_trk_p4 = mus_trk_p4_h.product();
  
  //muon type
  InputTag mus_type_tag(muonsInputTag.label(), "mustype");
  Handle<vector<int> > mus_type_h;
  iEvent.getByLabel(mus_type_tag, mus_type_h);
  const vector<int> *mus_type = mus_type_h.product();

  //-----------------------------------------------------------
  // electron variables
  //-----------------------------------------------------------
  InputTag els_charge_tag(electronsInputTag.label(),"elscharge");
  Handle<vector<int> > els_charge_h;
  iEvent.getByLabel(els_charge_tag, els_charge_h);
  const vector<int> *els_charge = els_charge_h.product();

  // electron p4
  InputTag els_p4_tag(electronsInputTag.label(),"elsp4");
  Handle<vector<LorentzVector> > els_p4_h;
  iEvent.getByLabel(els_p4_tag, els_p4_h);
  const vector<LorentzVector> *els_p4 = els_p4_h.product();
  
  //electrond0
  InputTag els_d0_tag(electronsInputTag.label(),"elsd0");
  Handle<vector<float> > els_d0_h;
  iEvent.getByLabel(els_d0_tag, els_d0_h);
  const vector<float> *els_d0 = els_d0_h.product();

  //electrond0 corrected from the beamSpot
  InputTag els_d0corr_tag(electronsInputTag.label(),"elsd0corr");
  Handle<vector<float> > els_d0corr_h;
  iEvent.getByLabel(els_d0corr_tag, els_d0corr_h);
  const vector<float> *els_d0corr = els_d0corr_h.product();

  //electron z0
  InputTag els_z0_tag(electronsInputTag.label(),"elsz0");
  Handle<vector<float> > els_z0_h;
  iEvent.getByLabel(els_z0_tag, els_z0_h);
  const vector<float> *els_z0 = els_z0_h.product();

  //electron z0, corrected for the beamspot
  InputTag els_z0corr_tag(electronsInputTag.label(),"elsz0corr");
  Handle<vector<float> > els_z0corr_h;
  iEvent.getByLabel(els_z0corr_tag, els_z0corr_h);
  const vector<float> *els_z0corr = els_z0corr_h.product();

  //chi2
  InputTag els_chi2_tag(electronsInputTag.label(),"elschi2");
  Handle<vector<float> > els_chi2_h;
  iEvent.getByLabel(els_chi2_tag, els_chi2_h);
  const vector<float> *els_chi2 = els_chi2_h.product();
  
  //ndof
  InputTag els_ndof_tag(electronsInputTag.label(),"elsndof");
  Handle<vector<float> > els_ndof_h;
  iEvent.getByLabel(els_ndof_tag, els_ndof_h);
  const vector<float> *els_ndof = els_ndof_h.product();
  
  //d0err
  InputTag els_d0err_tag(electronsInputTag.label(),"elsd0Err");
  Handle<vector<float> > els_d0err_h;
  iEvent.getByLabel(els_d0err_tag, els_d0err_h);
  const vector<float> *els_d0Err = els_d0err_h.product();

  //z0err
  InputTag els_z0err_tag(electronsInputTag.label(),"elsz0Err");
  Handle<vector<float> > els_z0err_h;
  iEvent.getByLabel(els_z0err_tag, els_z0err_h);
  const vector<float> *els_z0Err = els_z0err_h.product();

  //electron track P4
  InputTag els_trk_p4_tag(electronsInputTag.label(),"elstrkp4");
  Handle<vector<LorentzVector> > els_trk_p4_h;
  iEvent.getByLabel(els_trk_p4_tag, els_trk_p4_h);
  const vector<LorentzVector> *els_trk_p4 = els_trk_p4_h.product();

  //--------------------------------------------------------------------
  //Get the Jet collections
  //--------------------------------------------------------------------
  
  //jet p4
  Handle<vector<LorentzVector> > jets_p4_h;
  iEvent.getByLabel(jetsInputTag, jets_p4_h);
  const vector<LorentzVector> *jets_p4 = jets_p4_h.product();
 
  //event met
  InputTag met_tag(metInputTag.label(), "evtpfmet");
  Handle<float> met_tag_h;
  iEvent.getByLabel(met_tag, met_tag_h);
  const float* evt_met = met_tag_h.product();

  //event metPhi
  InputTag metphi_tag(metInputTag.label(), "evtpfmetPhi");
  Handle<float> metphi_tag_h;
  iEvent.getByLabel(metphi_tag, metphi_tag_h);
  const float* evt_metphi = metphi_tag_h.product();



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
            

      //fill the Jet vars
      vector<int> temp_jets_idx;
      vector<int> temp_other_jets_idx;
      vector<LorentzVector>  temp_jets_p4;      
      vector<LorentzVector>  temp_other_jets_p4;

      vector<LorentzVector> jets_nolep_p4;
	
      for(unsigned int i = 0; i<jets_p4->size(); i++) {
	
	// we don't want jets that overlap with electrons
	bool overlapsWithLepton = false;
	if(!testJetForLeptons(jets_p4->at(i), mus_p4->at(loose_index))) 
	  overlapsWithLepton = true;
	if(!testJetForLeptons(jets_p4->at(i), mus_p4->at(tight_index))) 
	  overlapsWithLepton = true;

	if( !overlapsWithLepton )
	  jets_nolep_p4        .push_back(jets_p4    ->at(i));
	
	double jet_eta = jets_p4->at(i).eta();
	double jet_pt  = jets_p4->at(i).Pt();
	
	if( fabs(jet_eta) < hypJetMaxEtaCut && jet_pt  > hypJetMinPtCut && !overlapsWithLepton) { //hyp jetas
	  temp_jets_idx.push_back(i);
	  temp_jets_p4                     .push_back(jets_p4              ->at(i));
	}
	else {
	  temp_other_jets_idx.push_back(i);
	  temp_other_jets_p4               .push_back(jets_p4              ->at(i));
	}
      }

      float temp_dPhi_nJet_Met = 9999.;
      float temp_sumJetPt = 0;

      for( unsigned int jidx = 0; jidx < temp_jets_p4.size(); jidx++ ) {

	float dphi       = fabs( temp_jets_p4[jidx].phi() - *evt_metphi          ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_metphi          : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_metphi          );

	if( dphi < temp_dPhi_nJet_Met )
	  temp_dPhi_nJet_Met = dphi;

	temp_sumJetPt += temp_jets_p4[jidx].pt();
      }
                                  
      float temp_Ht = temp_sumJetPt;

      temp_Ht += ( mus_p4->at(tight_index).pt() + mus_p4->at(loose_index).pt() + *evt_met );


      //push these into the hyp_jets and hyp_other_jets vars
      hyp_type          ->push_back(0                                     );
      hyp_p4            ->push_back(mus_p4->at(tight_index)+mus_p4->at(loose_index)               );
      hyp_lt_charge       ->push_back(mus_charge       ->at(tight_index)  );
      hyp_lt_index        ->push_back(tight_index                         );
      hyp_lt_id           ->push_back(-13*(mus_charge   ->at(tight_index)));
      hyp_lt_d0           ->push_back(mus_d0           ->at(tight_index)  );
      hyp_lt_z0           ->push_back(mus_z0           ->at(tight_index)  );
      hyp_lt_d0corr       ->push_back(mus_d0corr       ->at(tight_index)  );
      hyp_lt_z0corr       ->push_back(mus_z0corr       ->at(tight_index)  );
      hyp_lt_chi2         ->push_back(mus_chi2         ->at(tight_index)  );
      hyp_lt_ndof         ->push_back(mus_ndof         ->at(tight_index)  );
      hyp_lt_d0Err        ->push_back(mus_d0Err        ->at(tight_index)  );
      hyp_lt_z0Err        ->push_back(mus_z0Err        ->at(tight_index)  );
      hyp_lt_p4           ->push_back(mus_p4           ->at(tight_index)  );
      hyp_lt_trk_p4       ->push_back(mus_trk_p4       ->at(tight_index)  );
      hyp_ll_charge       ->push_back(mus_charge       ->at(loose_index)  );
      hyp_ll_index        ->push_back(loose_index                         );
      hyp_ll_id           ->push_back(-13*(mus_charge   ->at(loose_index)));
      hyp_ll_d0           ->push_back(mus_d0           ->at(loose_index)  );
      hyp_ll_z0           ->push_back(mus_z0           ->at(loose_index)  );
      hyp_ll_d0corr       ->push_back(mus_d0corr       ->at(loose_index)  );
      hyp_ll_z0corr       ->push_back(mus_z0corr       ->at(loose_index)  );
      hyp_ll_chi2         ->push_back(mus_chi2         ->at(loose_index)  );
      hyp_ll_ndof         ->push_back(mus_ndof         ->at(loose_index)  );
      hyp_ll_d0Err        ->push_back(mus_d0Err        ->at(loose_index)  );
      hyp_ll_z0Err        ->push_back(mus_z0Err        ->at(loose_index)  );
      hyp_ll_p4           ->push_back(mus_p4           ->at(loose_index)  );
      hyp_ll_trk_p4       ->push_back(mus_trk_p4       ->at(loose_index)  );
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
            
      //fill the Jet vars
      vector<int> temp_jets_idx;
      vector<int> temp_other_jets_idx;
      vector<LorentzVector>  temp_jets_p4;
      vector<LorentzVector>  temp_other_jets_p4;
      
      //these are for the MET correction later
      vector<LorentzVector> jets_nolep_p4;
      
      for(unsigned int i = 0; i<jets_p4->size(); i++) {

	// we don't want jets that overlap with electrons
	bool overlapsWithLepton = false;
	if(!testJetForLeptons(jets_p4->at(i), els_p4->at(loose_index))) 
	  overlapsWithLepton = true;
	if(!testJetForLeptons(jets_p4->at(i), els_p4->at(tight_index))) 
	  overlapsWithLepton = true;
	
	if( !overlapsWithLepton )
	  jets_nolep_p4        .push_back(jets_p4    ->at(i));
	
	double jet_eta = jets_p4->at(i).eta();
	double jet_pt = jets_p4->at(i).Pt();
	
	if( fabs(jet_eta) < hypJetMaxEtaCut && jet_pt  > hypJetMinPtCut && !overlapsWithLepton) { //hyp jets
	  temp_jets_idx.push_back(i);
	  temp_jets_p4                     .push_back(jets_p4              ->at(i));
	}
	else {
	  temp_other_jets_idx.push_back(i);
	  temp_other_jets_p4               .push_back(jets_p4              ->at(i));
	}//hyp_other Jets
	
      }//jet loop

      float temp_dPhi_nJet_Met = 9999.;

      float temp_sumJetPt = 0;

      for( unsigned int jidx = 0; jidx < temp_jets_p4.size(); jidx++ ) {

	float dphi       = fabs( temp_jets_p4[jidx].phi() - *evt_metphi            ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_metphi            : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_metphi            );

	if( dphi < temp_dPhi_nJet_Met )
	  temp_dPhi_nJet_Met = dphi;

	temp_sumJetPt += temp_jets_p4[jidx].pt();
      }

      float temp_Ht = temp_sumJetPt;

      temp_Ht += ( els_p4->at(tight_index).pt() + els_p4->at(loose_index).pt() + *evt_met );

      //hyp_Ht->push_back(temp_Ht);
    
      hyp_type          ->push_back(3);
      hyp_p4            ->push_back(els_p4->at(tight_index)+els_p4->at(loose_index)               );
      hyp_lt_charge       ->push_back(els_charge       ->at(tight_index)  );
      hyp_lt_index        ->push_back(tight_index                         );
      hyp_lt_id           ->push_back(-11*(els_charge   ->at(tight_index)));
      hyp_lt_d0           ->push_back(els_d0           ->at(tight_index)  );
      hyp_lt_z0           ->push_back(els_z0           ->at(tight_index)  );
      hyp_lt_d0corr       ->push_back(els_d0corr       ->at(tight_index)  );
      hyp_lt_z0corr       ->push_back(els_z0corr       ->at(tight_index)  );
      hyp_lt_chi2         ->push_back(els_chi2         ->at(tight_index)  );
      hyp_lt_ndof         ->push_back(els_ndof         ->at(tight_index)  );
      hyp_lt_d0Err        ->push_back(els_d0Err        ->at(tight_index)  );
      hyp_lt_z0Err        ->push_back(els_z0Err        ->at(tight_index)  );
      hyp_lt_p4           ->push_back(els_p4           ->at(tight_index)  );
      hyp_lt_trk_p4       ->push_back(els_trk_p4       ->at(tight_index)  );
      hyp_ll_charge       ->push_back(els_charge       ->at(loose_index)  );
      hyp_ll_index        ->push_back(loose_index                         );
      hyp_ll_id           ->push_back(-11*(els_charge   ->at(loose_index)));
      hyp_ll_d0           ->push_back(els_d0           ->at(loose_index)  );
      hyp_ll_z0           ->push_back(els_z0           ->at(loose_index)  );
      hyp_ll_d0corr       ->push_back(els_d0corr       ->at(loose_index)  );
      hyp_ll_z0corr       ->push_back(els_z0corr       ->at(loose_index)  );
      hyp_ll_chi2         ->push_back(els_chi2         ->at(loose_index)  );
      hyp_ll_ndof         ->push_back(els_ndof         ->at(loose_index)  );
      hyp_ll_d0Err        ->push_back(els_d0Err        ->at(loose_index)  );
      hyp_ll_z0Err        ->push_back(els_z0Err        ->at(loose_index)  );
      hyp_ll_p4           ->push_back(els_p4           ->at(loose_index)  );
      hyp_ll_trk_p4       ->push_back(els_trk_p4       ->at(loose_index)  );
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
                  
      //fill the Jet vars
      vector<int> temp_jets_idx;
      vector<int> temp_other_jets_idx;
      vector<LorentzVector>  temp_jets_p4;
      vector<LorentzVector>  temp_other_jets_p4;
      
      //these are for the MET correction later
      vector<LorentzVector> jets_nolep_p4;
       
      for(unsigned int i = 0; i<jets_p4->size(); i++) {
	
	bool overlapsWithLepton = false;
	//we don't any jets that overlap with an electron
	if(!testJetForLeptons(jets_p4->at(i), els_p4->at(els_index))) 
	  overlapsWithLepton = true;
	if(!testJetForLeptons(jets_p4->at(i), mus_p4->at(mus_index))) 
	  overlapsWithLepton = true;
	
	if( !overlapsWithLepton )
	  jets_nolep_p4        .push_back(jets_p4    ->at(i));

	double jet_eta = jets_p4->at(i).eta();
	double jet_pt  = jets_p4->at(i).Pt();

	
	if( fabs(jet_eta) < hypJetMaxEtaCut && jet_pt  > hypJetMinPtCut && !overlapsWithLepton) { 
	  temp_jets_idx.push_back(i);
	  temp_jets_p4                     .push_back(jets_p4              ->at(i));
	}
	else {	  
	  temp_other_jets_idx.push_back(i);
	  temp_other_jets_p4               .push_back(jets_p4              ->at(i));	  
	}
      }       

      hyp_p4            ->push_back(mus_p4->at(mus_index)+
				    els_p4->at(els_index)                 );
	
      if(el_pt < tightptcut && mu_pt > tightptcut) {
	hyp_type            ->push_back(1);
	  
	hyp_lt_charge       ->push_back(mus_charge       ->at(mus_index)  );
	hyp_lt_index        ->push_back(mus_index                         );
	hyp_lt_id           ->push_back(-13*(mus_charge   ->at(mus_index)));
	hyp_lt_d0           ->push_back(mus_d0           ->at(mus_index)  );
	hyp_lt_z0           ->push_back(mus_z0           ->at(mus_index)  );
	hyp_lt_d0corr       ->push_back(mus_d0corr       ->at(mus_index)  );
	hyp_lt_z0corr       ->push_back(mus_z0corr       ->at(mus_index)  );
	hyp_lt_chi2         ->push_back(mus_chi2         ->at(mus_index)  );
	hyp_lt_ndof         ->push_back(mus_ndof         ->at(mus_index)  );
	hyp_lt_d0Err        ->push_back(mus_d0Err        ->at(mus_index)  );
	hyp_lt_z0Err        ->push_back(mus_z0Err        ->at(mus_index)  );
	hyp_lt_p4           ->push_back(mus_p4           ->at(mus_index)  );
	hyp_lt_trk_p4       ->push_back(mus_trk_p4       ->at(mus_index)  );
	hyp_ll_charge       ->push_back(els_charge       ->at(els_index)  );
	hyp_ll_index        ->push_back(els_index                         );
	hyp_ll_id           ->push_back(-11*(els_charge   ->at(els_index)));
	hyp_ll_d0           ->push_back(els_d0           ->at(els_index)  );
	hyp_ll_z0           ->push_back(els_z0           ->at(els_index)  );
	hyp_ll_d0corr       ->push_back(els_d0corr       ->at(els_index)  );
	hyp_ll_z0corr       ->push_back(els_z0corr       ->at(els_index)  );
	hyp_ll_chi2         ->push_back(els_chi2         ->at(els_index)  );
	hyp_ll_ndof         ->push_back(els_ndof         ->at(els_index)  );
	hyp_ll_d0Err        ->push_back(els_d0Err        ->at(els_index)  );
	hyp_ll_z0Err        ->push_back(els_z0Err        ->at(els_index)  );
	hyp_ll_p4           ->push_back(els_p4           ->at(els_index)  );
	hyp_ll_trk_p4       ->push_back(els_trk_p4       ->at(els_index)  );
      }
    else {
	hyp_type            ->push_back(2);
	hyp_lt_charge       ->push_back(els_charge       ->at(els_index)  );
	hyp_lt_index        ->push_back(els_index                         );
	hyp_lt_id           ->push_back(-11*(els_charge   ->at(els_index)));
	hyp_lt_d0           ->push_back(els_d0           ->at(els_index)  );
	hyp_lt_z0           ->push_back(els_z0           ->at(els_index)  );
	hyp_lt_d0corr       ->push_back(els_d0corr       ->at(els_index)  );
	hyp_lt_z0corr       ->push_back(els_z0corr       ->at(els_index)  );
	hyp_lt_chi2         ->push_back(els_chi2         ->at(els_index)  );
	hyp_lt_ndof         ->push_back(els_ndof         ->at(els_index)  );
	hyp_lt_d0Err        ->push_back(els_d0Err        ->at(els_index)  );
	hyp_lt_z0Err        ->push_back(els_z0Err        ->at(els_index)  );
	hyp_lt_p4           ->push_back(els_p4           ->at(els_index)  );
	hyp_lt_trk_p4       ->push_back(els_trk_p4       ->at(els_index)  );
	
	hyp_ll_charge       ->push_back(mus_charge       ->at(mus_index)  );
	hyp_ll_index        ->push_back(mus_index                         );
	hyp_ll_id           ->push_back(-13*(mus_charge   ->at(mus_index)));
	hyp_ll_d0           ->push_back(mus_d0           ->at(mus_index)  );
	hyp_ll_z0           ->push_back(mus_z0           ->at(mus_index)  );
	hyp_ll_d0corr       ->push_back(mus_d0corr       ->at(mus_index)  );
	hyp_ll_z0corr       ->push_back(mus_z0corr       ->at(mus_index)  );
	hyp_ll_chi2         ->push_back(mus_chi2         ->at(mus_index)  );
	hyp_ll_ndof         ->push_back(mus_ndof         ->at(mus_index)  );
	hyp_ll_d0Err        ->push_back(mus_d0Err        ->at(mus_index)  );
	hyp_ll_z0Err        ->push_back(mus_z0Err        ->at(mus_index)  );
	hyp_ll_p4           ->push_back(mus_p4           ->at(mus_index)  );
	hyp_ll_trk_p4       ->push_back(mus_trk_p4       ->at(mus_index)  );
	
      }
    }
  }

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(hyp_type                     ,branchprefix+"type"                     );
  iEvent.put(hyp_p4                       ,branchprefix+"p4"                       );
  iEvent.put(hyp_lt_charge                ,branchprefix+"ltcharge"                 );
  iEvent.put(hyp_lt_index                 ,branchprefix+"ltindex"                  );
  iEvent.put(hyp_lt_id                    ,branchprefix+"ltid"                     );
  iEvent.put(hyp_lt_d0                    ,branchprefix+"ltd0"                     );
  iEvent.put(hyp_lt_z0                    ,branchprefix+"ltz0"                     );
  iEvent.put(hyp_lt_d0corr                ,branchprefix+"ltd0corr"                 );
  iEvent.put(hyp_lt_z0corr                ,branchprefix+"ltz0corr"                 );
  iEvent.put(hyp_lt_chi2                  ,branchprefix+"ltchi2"                   );
  iEvent.put(hyp_lt_ndof                  ,branchprefix+"ltndof"                   );
  iEvent.put(hyp_lt_d0Err                 ,branchprefix+"ltd0Err"                  );
  iEvent.put(hyp_lt_z0Err                 ,branchprefix+"ltz0Err"                  );
  iEvent.put(hyp_lt_p4                    ,branchprefix+"ltp4"                     );
  iEvent.put(hyp_lt_trk_p4                ,branchprefix+"lttrkp4"                  );
  iEvent.put(hyp_ll_charge                ,branchprefix+"llcharge"                 );
  iEvent.put(hyp_ll_index                 ,branchprefix+"llindex"                  );
  iEvent.put(hyp_ll_id                    ,branchprefix+"llid"                     );
  iEvent.put(hyp_ll_d0                    ,branchprefix+"lld0"                     );
  iEvent.put(hyp_ll_z0                    ,branchprefix+"llz0"                     );
  iEvent.put(hyp_ll_d0corr                ,branchprefix+"lld0corr"                 );
  iEvent.put(hyp_ll_z0corr                ,branchprefix+"llz0corr"                 );
  iEvent.put(hyp_ll_chi2                  ,branchprefix+"llchi2"                   );
  iEvent.put(hyp_ll_ndof                  ,branchprefix+"llndof"                   );
  iEvent.put(hyp_ll_d0Err                 ,branchprefix+"lld0Err"                  );
  iEvent.put(hyp_ll_z0Err                 ,branchprefix+"llz0Err"                  );
  iEvent.put(hyp_ll_p4                    ,branchprefix+"llp4"                     );
  iEvent.put(hyp_ll_trk_p4                ,branchprefix+"lltrkp4"                  );
}

// ------------ method called once each job just before starting event loop  ------------
void HypDilepMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HypDilepMaker::endJob() {
}

//----------------------------------------------------------------------------------------
bool HypDilepMaker::testJetForLeptons(const LorentzVector& jetP4, const LorentzVector& lepp4) {
  
  
  bool matched = false;
  float lepphi = lepp4.Phi();
  float jetphi = jetP4.Phi();
   
  float lepeta = lepp4.Eta();
  float jeteta = jetP4.Eta();
   
  float dphi = lepphi - jetphi;
  float deta = lepeta - jeteta;
  if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
   
  double dR = sqrt(dphi*dphi + deta*deta);
  if (dR < 0.4) 
    matched = true;
  
  return !matched;
}

//wrapper to use the mt2 class that spencer wrote. Could get rid of this entirely, and modify the mt2 class's functions but I'm too lazy.
double HypDilepMaker::mT2_bisect(const LorentzVector lep1_p4, const LorentzVector lep2_p4, 
				 const double met, const double metPhi) { 

  double pa[3];
  double pb[3];
  double pmiss[3];

  pa[0] = lep1_p4.M();
  pa[1] = lep1_p4.Px();
  pa[2] = lep1_p4.Py();

  pb[0] = lep2_p4.M();
  pb[1] = lep2_p4.Px();
  pb[2] = lep2_p4.Py();

  pmiss[0] = 0.;
  pmiss[1] = met*cos(metPhi);
  pmiss[2] = met*sin(metPhi);

  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta( pa, pb, pmiss );
  mt2_event.set_mn(0);

  double mt2_value = mt2_event.get_mt2();

  return mt2_value;
}


//define this as a plug-in
DEFINE_FWK_MODULE(HypDilepMaker);

