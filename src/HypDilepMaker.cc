// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      HypDilepMaker
// 
/**\class HypDilepMaker HypDilepMaker.cc CMS2/NtupleMaker/src/HypDilepMaker.cc

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
//
// Original Author:  Puneeth Kalavase
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: HypDilepMaker.cc,v 1.18 2009/09/02 20:05:37 fgolf Exp $
//
//

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
#include "CMS2/NtupleMaker/interface/HypDilepMaker.h"
#include "CMS2/NtupleMaker/interface/MatchUtilities.h"

#include "CMS2/NtupleMaker/interface/MT2Utility.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "Math/VectorUtil.h"
#include "TMath.h"


typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;


HypDilepMaker::HypDilepMaker(const edm::ParameterSet& iConfig)
{
  
  muonsInputTag            = iConfig.getParameter<InputTag>("muonsInputTag"                                        );
  electronsInputTag        = iConfig.getParameter<InputTag>("electronsInputTag"                                    );
  tcmetInputTag            = iConfig.getParameter<InputTag>("tcmetInputTag"                                        );
  metInputTag              = iConfig.getParameter<InputTag>("metInputTag"                                          );
  jetsInputTag             = iConfig.getParameter<InputTag>("jetsInputTag"                                         );
  trksInputTag             = iConfig.getParameter<InputTag>("trksInputTag"                                         );
  hypJetMaxEtaCut          = iConfig.getParameter<double>  ("hypJetMaxEtaCut"                                      );
  hypJetMinPtCut           = iConfig.getParameter<double>  ("hypJetMinPtCut"                                       );
  tightptcut               = iConfig.getParameter<double>  ("TightLepton_PtCut"                                    );
  looseptcut               = iConfig.getParameter<double>  ("LooseLepton_PtCut"                                    );

  produces<vector<int> >           ("hyptype"                    ).setBranchAlias("hyp_type"                       );
  produces<vector<int> >           ("hypnjets"                   ).setBranchAlias("hyp_njets"                      );
  produces<vector<int> >           ("hypnojets"                  ).setBranchAlias("hyp_nojets"                     );  
  produces<vector<LorentzVector> > ("hypp4"                      ).setBranchAlias("hyp_p4"                         );
  
  produces<vector<int> >           ("hypltvalidHits"             ).setBranchAlias("hyp_lt_validHits"               );
  produces<vector<int> >           ("hypltlostHits"              ).setBranchAlias("hyp_lt_lostHits"                );
  produces<vector<int> >           ("hypltcharge"                ).setBranchAlias("hyp_lt_charge"                  );
  produces<vector<int> >           ("hypltindex"                 ).setBranchAlias("hyp_lt_index"                   );
  produces<vector<int> >           ("hypltid"                    ).setBranchAlias("hyp_lt_id"                      );
  produces<vector<float> >         ("hypltd0"                    ).setBranchAlias("hyp_lt_d0"                      );
  produces<vector<float> >         ("hypltz0"                    ).setBranchAlias("hyp_lt_z0"                      );
  produces<vector<float> >         ("hypltd0corr"                ).setBranchAlias("hyp_lt_d0corr"                  );
  produces<vector<float> >         ("hypltz0corr"                ).setBranchAlias("hyp_lt_z0corr"                  );
  produces<vector<float> >         ("hypltchi2"                  ).setBranchAlias("hyp_lt_chi2"                    );
  produces<vector<float> >         ("hypltndof"                  ).setBranchAlias("hyp_lt_ndof"                    );
  produces<vector<float> >         ("hypltd0Err"                 ).setBranchAlias("hyp_lt_d0Err"                   );
  produces<vector<float> >         ("hypltz0Err"                 ).setBranchAlias("hyp_lt_z0Err"                   );
  produces<vector<float> >         ("hypltptErr"                 ).setBranchAlias("hyp_lt_ptErr"                   );
  produces<vector<float> >         ("hypltetaErr"                ).setBranchAlias("hyp_lt_etaErr"                  );
  produces<vector<float> >         ("hypltphiErr"                ).setBranchAlias("hyp_lt_phiErr"                  );
  produces<vector<LorentzVector > >("hypltp4"                    ).setBranchAlias("hyp_lt_p4"                      );
  produces<vector<LorentzVector > >("hyplttrkp4"                 ).setBranchAlias("hyp_lt_trk_p4"                  );
  
  produces<vector<int> >           ("hypllvalidHits"             ).setBranchAlias("hyp_ll_validHits"               );
  produces<vector<int> >           ("hyplllostHits"              ).setBranchAlias("hyp_ll_lostHits"                );
  produces<vector<int> >           ("hypllcharge"                ).setBranchAlias("hyp_ll_charge"                  );
  produces<vector<int> >           ("hypllindex"                 ).setBranchAlias("hyp_ll_index"                   );
  produces<vector<int> >           ("hypllid"                    ).setBranchAlias("hyp_ll_id"                      );
  produces<vector<float> >         ("hyplld0"                    ).setBranchAlias("hyp_ll_d0"                      );
  produces<vector<float> >         ("hypllz0"                    ).setBranchAlias("hyp_ll_z0"                      );
  produces<vector<float> >         ("hyplld0corr"                ).setBranchAlias("hyp_ll_d0corr"                  );
  produces<vector<float> >         ("hypllz0corr"                ).setBranchAlias("hyp_ll_z0corr"                  );
  produces<vector<float> >         ("hypllchi2"                  ).setBranchAlias("hyp_ll_chi2"                    );
  produces<vector<float> >         ("hypllndof"                  ).setBranchAlias("hyp_ll_ndof"                    );
  produces<vector<float> >         ("hyplld0Err"                 ).setBranchAlias("hyp_ll_d0Err"                   );
  produces<vector<float> >         ("hypllz0Err"                 ).setBranchAlias("hyp_ll_z0Err"                   );
  produces<vector<float> >         ("hypllptErr"                 ).setBranchAlias("hyp_ll_ptErr"                   );
  produces<vector<float> >         ("hyplletaErr"                ).setBranchAlias("hyp_ll_etaErr"                  );
  produces<vector<float> >         ("hypllphiErr"                ).setBranchAlias("hyp_ll_phiErr"                  );
  produces<vector<LorentzVector > >("hypllp4"                    ).setBranchAlias("hyp_ll_p4"                      );
  produces<vector<LorentzVector > >("hyplltrkp4"                 ).setBranchAlias("hyp_ll_trk_p4"                  );
  
  produces<vector<float> >         ("hypltdPhiunCorrMet"         ).setBranchAlias("hyp_lt_dPhi_unCorrMet"          );
  produces<vector<float> >         ("hyplldPhiunCorrMet"         ).setBranchAlias("hyp_ll_dPhi_unCorrMet"          );
  produces<vector<float> >         ("hypltdPhimuCorrMet"         ).setBranchAlias("hyp_lt_dPhi_muCorrMet"          );
  produces<vector<float> >         ("hyplldPhimuCorrMet"         ).setBranchAlias("hyp_ll_dPhi_muCorrMet"          );
  produces<vector<float> >         ("hypltdPhitcMet"             ).setBranchAlias("hyp_lt_dPhi_tcMet"              );
  produces<vector<float> >         ("hyplldPhitcMet"             ).setBranchAlias("hyp_ll_dPhi_tcMet"              );
  produces<vector<float> >         ("hypltdPhimetMuonJESCorr"    ).setBranchAlias("hyp_lt_dPhi_metMuonJESCorr"     );
  produces<vector<float> >         ("hyplldPhimetMuonJESCorr"    ).setBranchAlias("hyp_ll_dPhi_metMuonJESCorr"     );
  
  produces<vector<float> >         ("hypdPhinJetunCorrMet"       ).setBranchAlias("hyp_dPhi_nJet_unCorrMet"        );
  produces<vector<float> >         ("hypdPhinJetmuCorrMet"       ).setBranchAlias("hyp_dPhi_nJet_muCorrMet"        );
  produces<vector<float> >         ("hypdPhinJettcMet"           ).setBranchAlias("hyp_dPhi_nJet_tcMet"            );
  produces<vector<float> >         ("hypdPhinJetmetMuonJESCorr"  ).setBranchAlias("hyp_dPhi_nJet_metMuonJESCorr"   );
  
  produces<vector<float> >         ("hypsumJetPt"                ).setBranchAlias("hyp_sumJetPt"                   );
  produces<vector<float> >         ("hypHt"                      ).setBranchAlias("hyp_Ht"                         );
  
  //mt2
  produces<vector<float> >         ("hypmt2tcMet"                ).setBranchAlias("hyp_mt2_tcMet"                  );
  produces<vector<float> >         ("hypmt2muCorrMet"            ).setBranchAlias("hyp_mt2_muCorrMet"              );
  produces<vector<float> >         ("hypmt2metMuonJESCorr"       ).setBranchAlias("hyp_mt2_metMuonJESCorr"         );
  
  

  produces<vector<vector<int> > >  ("hypjetsidx"                 ).setBranchAlias("hyp_jets_idx"                   );
  produces<vector<vector<int> > >  ("hypotherjetsidx"            ).setBranchAlias("hyp_other_jets_idx"             );
  
  produces<vector<vector<LorentzVector> > >  ("hypjetsp4"       ).setBranchAlias("hyp_jets_p4"                     );
  produces<vector<vector<LorentzVector> > >  ("hypotherjetsp4"  ).setBranchAlias("hyp_other_jets_p4"               );
}


HypDilepMaker::~HypDilepMaker()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void HypDilepMaker::produce(Event& iEvent, const edm::EventSetup& iSetup)
{

  // output collections
  auto_ptr<vector<int> >           hyp_type                    (new vector<int>             );
  auto_ptr<vector<int> >           hyp_njets                   (new vector<int>             );
  auto_ptr<vector<int> >           hyp_nojets                  (new vector<int>             );
  auto_ptr<vector<LorentzVector> > hyp_p4                      (new vector<LorentzVector>   );

  auto_ptr<vector<int> >           hyp_lt_validHits            (new vector<int>             );
  auto_ptr<vector<int> >           hyp_lt_lostHits             (new vector<int>             );
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
  auto_ptr<vector<float> >         hyp_lt_ptErr                (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_etaErr               (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_phiErr               (new vector<float>           );
  auto_ptr<vector<LorentzVector> > hyp_lt_p4                   (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > hyp_lt_trk_p4               (new vector<LorentzVector>   );
  
  auto_ptr<vector<int> >           hyp_ll_validHits            (new vector<int>             );
  auto_ptr<vector<int> >           hyp_ll_lostHits             (new vector<int>             );
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
  auto_ptr<vector<float> >         hyp_ll_ptErr                (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_etaErr               (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_phiErr               (new vector<float>           );
  auto_ptr<vector<LorentzVector> > hyp_ll_p4                   (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > hyp_ll_trk_p4               (new vector<LorentzVector>   );
  
  auto_ptr<vector<float> >         hyp_lt_dPhi_unCorrMet       (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_dPhi_unCorrMet       (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_dPhi_muCorrMet       (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_dPhi_muCorrMet       (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_dPhi_tcMet           (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_dPhi_tcMet           (new vector<float>           );
  auto_ptr<vector<float> >         hyp_lt_dPhi_metMuonJESCorr  (new vector<float>           );
  auto_ptr<vector<float> >         hyp_ll_dPhi_metMuonJESCorr  (new vector<float>           );
  
  auto_ptr<vector<float> >         hyp_dPhi_nJet_unCorrMet     (new vector<float>           );
  auto_ptr<vector<float> >         hyp_dPhi_nJet_muCorrMet     (new vector<float>           );
  auto_ptr<vector<float> >         hyp_dPhi_nJet_tcMet         (new vector<float>           );
  auto_ptr<vector<float> >         hyp_dPhi_nJet_metMuonJESCorr(new vector<float>           );
  
  auto_ptr<vector<float> >         hyp_sumJetPt                (new vector<float>           );
  auto_ptr<vector<float> >         hyp_Ht                      (new vector<float>           );

  auto_ptr<vector<float> >         hyp_mt2_tcMet               (new vector<float>           );
  auto_ptr<vector<float> >         hyp_mt2_muCorrMet           (new vector<float>           );
  auto_ptr<vector<float> >         hyp_mt2_metMuonJESCorr      (new vector<float>           );
  

  auto_ptr<vector<vector<int> > >  hyp_jets_idx                (new vector<vector<int> >    );
  auto_ptr<vector<vector<int> > >  hyp_other_jets_idx          (new vector<vector<int> >    );
  auto_ptr<vector<vector<LorentzVector> > >  hyp_jets_p4       (new vector<vector<LorentzVector> > );
  auto_ptr<vector<vector<LorentzVector> > >  hyp_other_jets_p4 (new vector<vector<LorentzVector> > );
  
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

  //# of validHits on the muon track
  InputTag mus_validHits_tag(muonsInputTag.label(),"musvalidHits");
  Handle<vector<int> > mus_validHits_h;
  iEvent.getByLabel(mus_validHits_tag, mus_validHits_h);
  const vector<int> *mus_validHits = mus_validHits_h.product();

  //# of lostHits on the muon track
  InputTag mus_lostHits_tag(muonsInputTag.label(),"muslostHits");
  Handle<vector<int> > mus_lostHits_h;
  iEvent.getByLabel(mus_lostHits_tag, mus_lostHits_h);
  const vector<int> *mus_lostHits = mus_lostHits_h.product();
  
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

  //pterr
  InputTag mus_pterr_tag(muonsInputTag.label(),"musptErr");
  Handle<vector<float> > mus_pterr_h;
  iEvent.getByLabel(mus_pterr_tag, mus_pterr_h);
  const vector<float> *mus_ptErr = mus_pterr_h.product();

  //etaerr
  InputTag mus_etaerr_tag(muonsInputTag.label(),"musetaErr");
  Handle<vector<float> > mus_etaerr_h;
  iEvent.getByLabel(mus_etaerr_tag, mus_etaerr_h);
  const vector<float> *mus_etaErr = mus_etaerr_h.product();

  //phierr
  InputTag mus_phierr_tag(muonsInputTag.label(),"musphiErr");
  Handle<vector<float> > mus_phierr_h;
  iEvent.getByLabel(mus_phierr_tag, mus_phierr_h);
  const vector<float> *mus_phiErr = mus_phierr_h.product();

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
  
  //# of validHits on the el track
  InputTag els_validHits_tag(electronsInputTag.label(),"elsvalidHits");
  Handle<vector<int> > els_validHits_h;
  iEvent.getByLabel(els_validHits_tag, els_validHits_h);
  const vector<int> *els_validHits = els_validHits_h.product();

  //# of lostHits on the electron track
  InputTag els_lostHits_tag(electronsInputTag.label(),"elslostHits");
  Handle<vector<int> > els_lostHits_h;
  iEvent.getByLabel(els_lostHits_tag, els_lostHits_h);
  const vector<int> *els_lostHits = els_lostHits_h.product();

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

  //pterr
  InputTag els_pterr_tag(electronsInputTag.label(),"elsptErr");
  Handle<vector<float> > els_pterr_h;
  iEvent.getByLabel(els_pterr_tag, els_pterr_h);
  const vector<float> *els_ptErr = els_pterr_h.product();

  //etaerr
  InputTag els_etaerr_tag(electronsInputTag.label(),"elsetaErr");
  Handle<vector<float> > els_etaerr_h;
  iEvent.getByLabel(els_etaerr_tag, els_etaerr_h);
  const vector<float> *els_etaErr = els_etaerr_h.product();

  //phierr
  InputTag els_phierr_tag(electronsInputTag.label(),"elsphiErr");
  Handle<vector<float> > els_phierr_h;
  iEvent.getByLabel(els_phierr_tag, els_phierr_h);
  const vector<float> *els_phiErr = els_phierr_h.product();

  //electron track P4
  InputTag els_trk_p4_tag(electronsInputTag.label(),"elstrkp4");
  Handle<vector<LorentzVector> > els_trk_p4_h;
  iEvent.getByLabel(els_trk_p4_tag, els_trk_p4_h);
  const vector<LorentzVector> *els_trk_p4 = els_trk_p4_h.product();

  //--------------------------------------------------------------------
  //Get the Jet collections
  //--------------------------------------------------------------------
  
  //jet p4
  InputTag jets_p4_tag(jetsInputTag.label(), "jptsp4");
  Handle<vector<LorentzVector> > jets_p4_h;
  iEvent.getByLabel(jets_p4_tag, jets_p4_h);
  const vector<LorentzVector> *jets_p4 = jets_p4_h.product();
 
  //event metPhi
  InputTag metphi_tag(metInputTag.label(), "evtmetPhi");
  Handle<float> metphi_tag_h;
  iEvent.getByLabel(metphi_tag, metphi_tag_h);
  const float* evt_metphi = metphi_tag_h.product();

  //muon corrected met
  InputTag mumet_tag(metInputTag.label(), "evtmetMuonCorr");
  Handle<float> mumet_tag_h;
  iEvent.getByLabel(mumet_tag, mumet_tag_h);
  const float* evt_mumet = mumet_tag_h.product();

  //muon corrected metPhi
  InputTag mumetphi_tag(metInputTag.label(), "evtmetMuonCorrPhi");
  Handle<float> mumetphi_tag_h;
  iEvent.getByLabel(mumetphi_tag, mumetphi_tag_h);
  const float* evt_mumetphi = mumetphi_tag_h.product();

  //muon+JES corrected met
  InputTag muonJESCorrmet_tag(metInputTag.label(), "evtmetMuonJESCorr");
  Handle<float> muonJESCorrmet_tag_h;
  iEvent.getByLabel(muonJESCorrmet_tag, muonJESCorrmet_tag_h);
  const float* evt_metMuonJESCorr = muonJESCorrmet_tag_h.product();

  //muon+JES corrected metPhi
  InputTag muonJESCorrmetphi_tag(metInputTag.label(), "evtmetMuonJESCorrPhi");
  Handle<float> muonJESCorrmetphi_tag_h;
  iEvent.getByLabel(muonJESCorrmetphi_tag, muonJESCorrmetphi_tag_h);
  const float* evt_metMuonJESCorrPhi = muonJESCorrmetphi_tag_h.product();

  //event met
  InputTag tcmet_tag(tcmetInputTag.label(), "evttcmet");
  Handle<float> tcmet_tag_h;
  iEvent.getByLabel(tcmet_tag, tcmet_tag_h);
  const float* evt_tcmet = tcmet_tag_h.product();

  //event metPhi
  InputTag tcmetphi_tag(tcmetInputTag.label(), "evttcmetPhi");
  Handle<float> tcmetphi_tag_h;
  iEvent.getByLabel(tcmetphi_tag, tcmetphi_tag_h);
  const float* evt_tcmetphi = tcmetphi_tag_h.product();

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
      
      //mt2 shit, testing
      hyp_mt2_tcMet          ->push_back(mT2_bisect(mus_p4->at(tight_index), mus_p4->at(loose_index), 
						    *evt_tcmet, *evt_tcmetphi) );
      hyp_mt2_muCorrMet      ->push_back(mT2_bisect(mus_p4->at(tight_index), mus_p4->at(loose_index), 
						    *evt_mumet, *evt_mumetphi) );
      hyp_mt2_metMuonJESCorr ->push_back(mT2_bisect(mus_p4->at(tight_index), mus_p4->at(loose_index), 
						    *evt_metMuonJESCorr, *evt_metMuonJESCorrPhi ) );
      

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

      hyp_jets_idx->push_back(temp_jets_idx);
      hyp_other_jets_idx->push_back(temp_other_jets_idx);
      hyp_jets_p4                    ->push_back(temp_jets_p4                  );
      hyp_other_jets_p4              ->push_back(temp_other_jets_p4            );

      float temp_lt_dphi_unCorrMet;
      float temp_ll_dphi_unCorrMet;
      float temp_lt_dphi_muCorrMet;
      float temp_ll_dphi_muCorrMet;
      float temp_lt_dphi_metMuonJESCorr;
      float temp_ll_dphi_metMuonJESCorr;
      float temp_lt_dphi_tcMet;
      float temp_ll_dphi_tcMet;

      temp_lt_dphi_unCorrMet = fabs( mus_p4->at(tight_index).phi() - *evt_metphi          ) < TMath::Pi() ? mus_p4->at(tight_index).phi() - *evt_metphi          : 2*TMath::Pi() - ( mus_p4->at(tight_index).phi() - *evt_metphi          );
      temp_ll_dphi_unCorrMet = fabs( mus_p4->at(loose_index).phi() - *evt_metphi          ) < TMath::Pi() ? mus_p4->at(loose_index).phi() - *evt_metphi          : 2*TMath::Pi() - ( mus_p4->at(loose_index).phi() - *evt_metphi          );
      temp_lt_dphi_muCorrMet = fabs( mus_p4->at(tight_index).phi() - *evt_mumetphi  ) < TMath::Pi() ? mus_p4->at(tight_index).phi() - *evt_mumetphi  : 2*TMath::Pi() - ( mus_p4->at(tight_index).phi() - *evt_mumetphi  );
      temp_ll_dphi_muCorrMet = fabs( mus_p4->at(loose_index).phi() - *evt_mumetphi  ) < TMath::Pi() ? mus_p4->at(loose_index).phi() - *evt_mumetphi  : 2*TMath::Pi() - ( mus_p4->at(loose_index).phi() - *evt_mumetphi  );
      temp_lt_dphi_metMuonJESCorr  = fabs( mus_p4->at(tight_index).phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? mus_p4->at(tight_index).phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( mus_p4->at(tight_index).phi() - *evt_metMuonJESCorrPhi );
      temp_ll_dphi_metMuonJESCorr  = fabs( mus_p4->at(loose_index).phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? mus_p4->at(loose_index).phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( mus_p4->at(loose_index).phi() - *evt_metMuonJESCorrPhi );
      temp_lt_dphi_tcMet     = fabs( mus_p4->at(tight_index).phi() - *evt_tcmetphi        ) < TMath::Pi() ? mus_p4->at(tight_index).phi() - *evt_tcmetphi        : 2*TMath::Pi() - ( mus_p4->at(tight_index).phi() - *evt_tcmetphi        );
      temp_ll_dphi_tcMet     = fabs( mus_p4->at(loose_index).phi() - *evt_tcmetphi        ) < TMath::Pi() ? mus_p4->at(loose_index).phi() - *evt_tcmetphi        : 2*TMath::Pi() - ( mus_p4->at(loose_index).phi() - *evt_tcmetphi        );

      hyp_lt_dPhi_unCorrMet->push_back(temp_lt_dphi_unCorrMet); 
      hyp_ll_dPhi_unCorrMet->push_back(temp_ll_dphi_unCorrMet); 
      hyp_lt_dPhi_muCorrMet->push_back(temp_lt_dphi_muCorrMet); 
      hyp_ll_dPhi_muCorrMet->push_back(temp_ll_dphi_muCorrMet); 
      hyp_lt_dPhi_metMuonJESCorr ->push_back(temp_lt_dphi_metMuonJESCorr );  
      hyp_ll_dPhi_metMuonJESCorr ->push_back(temp_ll_dphi_metMuonJESCorr );  
      hyp_lt_dPhi_tcMet    ->push_back(temp_lt_dphi_tcMet    );    
      hyp_ll_dPhi_tcMet    ->push_back(temp_ll_dphi_tcMet    );     

      float temp_dPhi_nJet_unCorrMet = 999;
      float temp_dPhi_nJet_muCorrMet = 999;
      float temp_dPhi_nJet_tcMet     = 999;   
      float temp_dPhi_nJet_metMuonJESCorr  = 999;

      float temp_sumJetPt = 0;

      for( unsigned int jidx = 0; jidx < temp_jets_p4.size(); jidx++ ) {

	float dphi       = fabs( temp_jets_p4[jidx].phi() - *evt_metphi          ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_metphi          : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_metphi          );
	float dphi_mu    = fabs( temp_jets_p4[jidx].phi() - *evt_mumetphi  ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_mumetphi  : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_mumetphi  );
	float dphi_muonJESCorr = fabs( temp_jets_p4[jidx].phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_metMuonJESCorrPhi );
	float dphi_tc    = fabs( temp_jets_p4[jidx].phi() - *evt_tcmetphi        ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_tcmetphi        : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_tcmetphi        );

	if( dphi < temp_dPhi_nJet_unCorrMet )
	  temp_dPhi_nJet_unCorrMet = dphi;

	if( dphi_mu < temp_dPhi_nJet_muCorrMet )
	  temp_dPhi_nJet_muCorrMet = dphi_mu;

	if( dphi_muonJESCorr < temp_dPhi_nJet_metMuonJESCorr )
	  temp_dPhi_nJet_metMuonJESCorr = dphi_muonJESCorr;

	if( dphi_tc < temp_dPhi_nJet_tcMet )
	  temp_dPhi_nJet_tcMet = dphi_tc;

	temp_sumJetPt += temp_jets_p4[jidx].pt();
      }
                                  
      hyp_dPhi_nJet_unCorrMet->push_back(temp_dPhi_nJet_unCorrMet);
      hyp_dPhi_nJet_muCorrMet->push_back(temp_dPhi_nJet_muCorrMet);
      hyp_dPhi_nJet_tcMet    ->push_back(temp_dPhi_nJet_tcMet    );   
      hyp_dPhi_nJet_metMuonJESCorr ->push_back(temp_dPhi_nJet_metMuonJESCorr );

      hyp_sumJetPt->push_back( temp_sumJetPt );

      float temp_Ht = temp_sumJetPt;

      temp_Ht += ( mus_p4->at(tight_index).pt() + mus_p4->at(loose_index).pt() + *evt_tcmet );

      hyp_Ht->push_back(temp_Ht);

      //push these into the hyp_jets and hyp_other_jets vars
      hyp_type          ->push_back(0                                     );
      hyp_njets         ->push_back(temp_jets_p4.size()                   );     
      hyp_nojets        ->push_back(temp_other_jets_p4.size()             );     
      hyp_p4            ->push_back(mus_p4->at(tight_index)+mus_p4->at(loose_index)               );
    
      hyp_lt_validHits    ->push_back(mus_validHits    ->at(tight_index)  );
      hyp_lt_lostHits     ->push_back(mus_lostHits     ->at(tight_index)  );
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
      hyp_lt_ptErr        ->push_back(mus_ptErr        ->at(tight_index)  );
      hyp_lt_etaErr       ->push_back(mus_etaErr       ->at(tight_index)  );
      hyp_lt_phiErr       ->push_back(mus_phiErr       ->at(tight_index)  );
      hyp_lt_p4           ->push_back(mus_p4           ->at(tight_index)  );
      hyp_lt_trk_p4       ->push_back(mus_trk_p4       ->at(tight_index)  );
      
      hyp_ll_validHits    ->push_back(mus_validHits    ->at(loose_index)  );
      hyp_ll_lostHits     ->push_back(mus_lostHits     ->at(loose_index)  );
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
      hyp_ll_ptErr        ->push_back(mus_ptErr        ->at(loose_index)  );
      hyp_ll_etaErr       ->push_back(mus_etaErr       ->at(loose_index)  );
      hyp_ll_phiErr       ->push_back(mus_phiErr       ->at(loose_index)  );
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
      
      //mt2 shit, testing
      hyp_mt2_tcMet          ->push_back(mT2_bisect(els_p4->at(tight_index), els_p4->at(loose_index), 
						    *evt_tcmet, *evt_tcmetphi) );
      hyp_mt2_muCorrMet      ->push_back(mT2_bisect(els_p4->at(tight_index), els_p4->at(loose_index), 
						    *evt_mumet, *evt_mumetphi) );
      hyp_mt2_metMuonJESCorr ->push_back(mT2_bisect(els_p4->at(tight_index), els_p4->at(loose_index), 
						    *evt_metMuonJESCorr, *evt_metMuonJESCorrPhi ) );

      
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

      //push these into the hyp_jets and hyp_other_jets vars
      hyp_jets_idx->push_back(temp_jets_idx);
      hyp_other_jets_idx->push_back(temp_other_jets_idx);
      hyp_jets_p4                    ->push_back(temp_jets_p4                  );
      hyp_other_jets_p4              ->push_back(temp_other_jets_p4            );

      float temp_lt_dphi_unCorrMet;
      float temp_ll_dphi_unCorrMet;
      float temp_lt_dphi_muCorrMet;
      float temp_ll_dphi_muCorrMet;
      float temp_lt_dphi_metMuonJESCorr;
      float temp_ll_dphi_metMuonJESCorr;
      float temp_lt_dphi_tcMet;
      float temp_ll_dphi_tcMet;

      temp_lt_dphi_unCorrMet = fabs( els_p4->at(tight_index).phi() - *evt_metphi            ) < TMath::Pi() ? els_p4->at(tight_index).phi() - *evt_metphi            : 2*TMath::Pi() - ( els_p4->at(tight_index).phi() - *evt_metphi            );
      temp_ll_dphi_unCorrMet = fabs( els_p4->at(loose_index).phi() - *evt_metphi            ) < TMath::Pi() ? els_p4->at(loose_index).phi() - *evt_metphi            : 2*TMath::Pi() - ( els_p4->at(loose_index).phi() - *evt_metphi            );
      temp_lt_dphi_muCorrMet = fabs( els_p4->at(tight_index).phi() - *evt_mumetphi          ) < TMath::Pi() ? els_p4->at(tight_index).phi() - *evt_mumetphi          : 2*TMath::Pi() - ( els_p4->at(tight_index).phi() - *evt_mumetphi          );
      temp_ll_dphi_muCorrMet = fabs( els_p4->at(loose_index).phi() - *evt_mumetphi          ) < TMath::Pi() ? els_p4->at(loose_index).phi() - *evt_mumetphi          : 2*TMath::Pi() - ( els_p4->at(loose_index).phi() - *evt_mumetphi          );
      temp_lt_dphi_metMuonJESCorr  = fabs( els_p4->at(tight_index).phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? els_p4->at(tight_index).phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( els_p4->at(tight_index).phi() - *evt_metMuonJESCorrPhi );
      temp_ll_dphi_metMuonJESCorr  = fabs( els_p4->at(loose_index).phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? els_p4->at(loose_index).phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( els_p4->at(loose_index).phi() - *evt_metMuonJESCorrPhi );
      temp_lt_dphi_tcMet     = fabs( els_p4->at(tight_index).phi() - *evt_tcmetphi          ) < TMath::Pi() ? els_p4->at(tight_index).phi() - *evt_tcmetphi          : 2*TMath::Pi() - ( els_p4->at(tight_index).phi() - *evt_tcmetphi          );
      temp_ll_dphi_tcMet     = fabs( els_p4->at(loose_index).phi() - *evt_tcmetphi          ) < TMath::Pi() ? els_p4->at(loose_index).phi() - *evt_tcmetphi          : 2*TMath::Pi() - ( els_p4->at(loose_index).phi() - *evt_tcmetphi          );

      hyp_lt_dPhi_unCorrMet->push_back(temp_lt_dphi_unCorrMet); 
      hyp_ll_dPhi_unCorrMet->push_back(temp_ll_dphi_unCorrMet); 
      hyp_lt_dPhi_muCorrMet->push_back(temp_lt_dphi_muCorrMet); 
      hyp_ll_dPhi_muCorrMet->push_back(temp_ll_dphi_muCorrMet); 
      hyp_lt_dPhi_metMuonJESCorr ->push_back(temp_lt_dphi_metMuonJESCorr );  
      hyp_ll_dPhi_metMuonJESCorr ->push_back(temp_ll_dphi_metMuonJESCorr );  
      hyp_lt_dPhi_tcMet    ->push_back(temp_lt_dphi_tcMet    );    
      hyp_ll_dPhi_tcMet    ->push_back(temp_ll_dphi_tcMet    );     

      float temp_dPhi_nJet_unCorrMet = 999;
      float temp_dPhi_nJet_muCorrMet = 999;
      float temp_dPhi_nJet_tcMet     = 999;   
      float temp_dPhi_nJet_metMuonJESCorr  = 999;

      float temp_sumJetPt = 0;

      for( unsigned int jidx = 0; jidx < temp_jets_p4.size(); jidx++ ) {

	float dphi       = fabs( temp_jets_p4[jidx].phi() - *evt_metphi            ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_metphi            : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_metphi            );
	float dphi_mu    = fabs( temp_jets_p4[jidx].phi() - *evt_mumetphi          ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_mumetphi          : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_mumetphi          );
	float dphi_muonJESCorr = fabs( temp_jets_p4[jidx].phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_metMuonJESCorrPhi );
	float dphi_tc    = fabs( temp_jets_p4[jidx].phi() - *evt_tcmetphi          ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_tcmetphi          : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_tcmetphi          );

	if( dphi < temp_dPhi_nJet_unCorrMet )
	  temp_dPhi_nJet_unCorrMet = dphi;

	if( dphi_mu < temp_dPhi_nJet_muCorrMet )
	  temp_dPhi_nJet_muCorrMet = dphi_mu;

	if( dphi_muonJESCorr < temp_dPhi_nJet_metMuonJESCorr )
	  temp_dPhi_nJet_metMuonJESCorr = dphi_muonJESCorr;

	if( dphi_tc < temp_dPhi_nJet_tcMet )
	  temp_dPhi_nJet_tcMet = dphi_tc;

	temp_sumJetPt += temp_jets_p4[jidx].pt();
      }
                                  
      hyp_dPhi_nJet_unCorrMet->push_back(temp_dPhi_nJet_unCorrMet);
      hyp_dPhi_nJet_muCorrMet->push_back(temp_dPhi_nJet_muCorrMet);
      hyp_dPhi_nJet_tcMet    ->push_back(temp_dPhi_nJet_tcMet    );   
      hyp_dPhi_nJet_metMuonJESCorr ->push_back(temp_dPhi_nJet_metMuonJESCorr );

      hyp_sumJetPt->push_back( temp_sumJetPt );

      float temp_Ht = temp_sumJetPt;

      temp_Ht += ( els_p4->at(tight_index).pt() + els_p4->at(loose_index).pt() + *evt_tcmet );

      hyp_Ht->push_back(temp_Ht);
    
      hyp_type          ->push_back(3);
      hyp_njets         ->push_back(temp_jets_p4.size()                   );     
      hyp_nojets        ->push_back(temp_other_jets_p4.size()             );     
      hyp_p4            ->push_back(els_p4->at(tight_index)+
				    els_p4->at(loose_index)               );
      
      hyp_lt_validHits    ->push_back(els_validHits    ->at(tight_index)  );
      hyp_lt_lostHits     ->push_back(els_lostHits     ->at(tight_index)  );
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
      hyp_lt_ptErr        ->push_back(els_ptErr        ->at(tight_index)  );
      hyp_lt_etaErr       ->push_back(els_etaErr       ->at(tight_index)  );
      hyp_lt_phiErr       ->push_back(els_phiErr       ->at(tight_index)  );
      hyp_lt_p4           ->push_back(els_p4           ->at(tight_index)  );
      hyp_lt_trk_p4       ->push_back(els_trk_p4       ->at(tight_index)  );
      
      hyp_ll_validHits    ->push_back(els_validHits    ->at(loose_index)  );
      hyp_ll_lostHits     ->push_back(els_lostHits     ->at(loose_index)  );
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
      hyp_ll_ptErr        ->push_back(els_ptErr        ->at(loose_index)  );
      hyp_ll_etaErr       ->push_back(els_etaErr       ->at(loose_index)  );
      hyp_ll_phiErr       ->push_back(els_phiErr       ->at(loose_index)  );
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
      
      //mu is tight and e is loose only if el_pt < tightcut and 
      //mupt > tightcut      
      	
      //correct the met for the hyp muons
      pair<LorentzVector, LorentzVector> muon_pair = make_pair(mus_p4->at(mus_index),
							       mus_trk_p4->at(mus_index) );
      
      //mt2 shit, testing
      hyp_mt2_tcMet          ->push_back(mT2_bisect(els_p4->at(els_index), mus_p4->at(mus_index), 
						    *evt_tcmet, *evt_tcmetphi) );
      hyp_mt2_muCorrMet      ->push_back(mT2_bisect(els_p4->at(els_index), mus_p4->at(mus_index), 
						    *evt_mumet, *evt_mumetphi) );
      hyp_mt2_metMuonJESCorr ->push_back(mT2_bisect(els_p4->at(els_index), mus_p4->at(mus_index), 
						    *evt_metMuonJESCorr, *evt_metMuonJESCorrPhi ) );
      

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

      //push these into the hyp_jets and hyp_other_jets vars
      hyp_jets_idx->push_back(temp_jets_idx);
      hyp_other_jets_idx->push_back(temp_other_jets_idx);
      hyp_jets_p4                    ->push_back(temp_jets_p4                  );
      hyp_other_jets_p4              ->push_back(temp_other_jets_p4            );

      float temp_lt_dphi_unCorrMet;
      float temp_ll_dphi_unCorrMet;
      float temp_lt_dphi_muCorrMet;
      float temp_ll_dphi_muCorrMet;
      float temp_lt_dphi_metMuonJESCorr;
      float temp_ll_dphi_metMuonJESCorr;
      float temp_lt_dphi_tcMet;
      float temp_ll_dphi_tcMet;

      if(el_pt < tightptcut && mu_pt > tightptcut) {

	temp_lt_dphi_unCorrMet = fabs( mus_p4->at(mus_index).phi() - *evt_metphi            ) < TMath::Pi() ? mus_p4->at(mus_index).phi() - *evt_metphi            : 2*TMath::Pi() - ( mus_p4->at(mus_index).phi() - *evt_metphi            );
	temp_ll_dphi_unCorrMet = fabs( els_p4->at(els_index).phi() - *evt_metphi            ) < TMath::Pi() ? els_p4->at(els_index).phi() - *evt_metphi            : 2*TMath::Pi() - ( els_p4->at(els_index).phi() - *evt_metphi            );
	temp_lt_dphi_muCorrMet = fabs( mus_p4->at(mus_index).phi() - *evt_mumetphi          ) < TMath::Pi() ? mus_p4->at(mus_index).phi() - *evt_mumetphi          : 2*TMath::Pi() - ( mus_p4->at(mus_index).phi() - *evt_mumetphi          );
	temp_ll_dphi_muCorrMet = fabs( els_p4->at(els_index).phi() - *evt_mumetphi          ) < TMath::Pi() ? els_p4->at(els_index).phi() - *evt_mumetphi          : 2*TMath::Pi() - ( els_p4->at(els_index).phi() - *evt_mumetphi          );
	temp_lt_dphi_metMuonJESCorr  = fabs( mus_p4->at(mus_index).phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? mus_p4->at(mus_index).phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( mus_p4->at(mus_index).phi() - *evt_metMuonJESCorrPhi );
	temp_ll_dphi_metMuonJESCorr  = fabs( els_p4->at(els_index).phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? els_p4->at(els_index).phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( els_p4->at(els_index).phi() - *evt_metMuonJESCorrPhi );
	temp_lt_dphi_tcMet     = fabs( mus_p4->at(mus_index).phi() - *evt_tcmetphi          ) < TMath::Pi() ? mus_p4->at(mus_index).phi() - *evt_tcmetphi          : 2*TMath::Pi() - ( mus_p4->at(mus_index).phi() - *evt_tcmetphi          );
	temp_ll_dphi_tcMet     = fabs( els_p4->at(els_index).phi() - *evt_tcmetphi          ) < TMath::Pi() ? els_p4->at(els_index).phi() - *evt_tcmetphi          : 2*TMath::Pi() - ( els_p4->at(els_index).phi() - *evt_tcmetphi          );

	hyp_lt_dPhi_unCorrMet->push_back(temp_lt_dphi_unCorrMet); 
	hyp_ll_dPhi_unCorrMet->push_back(temp_ll_dphi_unCorrMet); 
	hyp_lt_dPhi_muCorrMet->push_back(temp_lt_dphi_muCorrMet); 
	hyp_ll_dPhi_muCorrMet->push_back(temp_ll_dphi_muCorrMet); 
	hyp_lt_dPhi_metMuonJESCorr ->push_back(temp_lt_dphi_metMuonJESCorr );  
	hyp_ll_dPhi_metMuonJESCorr ->push_back(temp_ll_dphi_metMuonJESCorr );  
	hyp_lt_dPhi_tcMet    ->push_back(temp_lt_dphi_tcMet    );    
	hyp_ll_dPhi_tcMet    ->push_back(temp_ll_dphi_tcMet    );     
      }

      else {

	temp_ll_dphi_unCorrMet = fabs( mus_p4->at(mus_index).phi() - *evt_metphi            ) < TMath::Pi() ? mus_p4->at(mus_index).phi() - *evt_metphi            : 2*TMath::Pi() - ( mus_p4->at(mus_index).phi() - *evt_metphi            );
	temp_lt_dphi_unCorrMet = fabs( els_p4->at(els_index).phi() - *evt_metphi            ) < TMath::Pi() ? els_p4->at(els_index).phi() - *evt_metphi            : 2*TMath::Pi() - ( els_p4->at(els_index).phi() - *evt_metphi            );
	temp_ll_dphi_muCorrMet = fabs( mus_p4->at(mus_index).phi() - *evt_mumetphi          ) < TMath::Pi() ? mus_p4->at(mus_index).phi() - *evt_mumetphi          : 2*TMath::Pi() - ( mus_p4->at(mus_index).phi() - *evt_mumetphi          );
	temp_lt_dphi_muCorrMet = fabs( els_p4->at(els_index).phi() - *evt_mumetphi          ) < TMath::Pi() ? els_p4->at(els_index).phi() - *evt_mumetphi          : 2*TMath::Pi() - ( els_p4->at(els_index).phi() - *evt_mumetphi          );
	temp_ll_dphi_metMuonJESCorr  = fabs( mus_p4->at(mus_index).phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? mus_p4->at(mus_index).phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( mus_p4->at(mus_index).phi() - *evt_metMuonJESCorrPhi );
	temp_lt_dphi_metMuonJESCorr  = fabs( els_p4->at(els_index).phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? els_p4->at(els_index).phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( els_p4->at(els_index).phi() - *evt_metMuonJESCorrPhi );
	temp_ll_dphi_tcMet     = fabs( mus_p4->at(mus_index).phi() - *evt_tcmetphi          ) < TMath::Pi() ? mus_p4->at(mus_index).phi() - *evt_tcmetphi          : 2*TMath::Pi() - ( mus_p4->at(mus_index).phi() - *evt_tcmetphi          );
	temp_lt_dphi_tcMet     = fabs( els_p4->at(els_index).phi() - *evt_tcmetphi          ) < TMath::Pi() ? els_p4->at(els_index).phi() - *evt_tcmetphi          : 2*TMath::Pi() - ( els_p4->at(els_index).phi() - *evt_tcmetphi          );

	hyp_lt_dPhi_unCorrMet->push_back(temp_lt_dphi_unCorrMet); 
	hyp_ll_dPhi_unCorrMet->push_back(temp_ll_dphi_unCorrMet); 
	hyp_lt_dPhi_muCorrMet->push_back(temp_lt_dphi_muCorrMet); 
	hyp_ll_dPhi_muCorrMet->push_back(temp_ll_dphi_muCorrMet); 
	hyp_lt_dPhi_metMuonJESCorr ->push_back(temp_lt_dphi_metMuonJESCorr );  
	hyp_ll_dPhi_metMuonJESCorr ->push_back(temp_ll_dphi_metMuonJESCorr );  
	hyp_lt_dPhi_tcMet    ->push_back(temp_lt_dphi_tcMet    );    
	hyp_ll_dPhi_tcMet    ->push_back(temp_ll_dphi_tcMet    );     
      }

      float temp_dPhi_nJet_unCorrMet = 999;
      float temp_dPhi_nJet_muCorrMet = 999;
      float temp_dPhi_nJet_tcMet     = 999;   
      float temp_dPhi_nJet_metMuonJESCorr  = 999;

      float temp_sumJetPt = 0;

      for( unsigned int jidx = 0; jidx < temp_jets_p4.size(); jidx++ ) {

	float dphi       = fabs( temp_jets_p4[jidx].phi() - *evt_metphi            ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_metphi            : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_metphi            );
	float dphi_mu    = fabs( temp_jets_p4[jidx].phi() - *evt_mumetphi          ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_mumetphi          : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_mumetphi          );
	float dphi_muonJESCorr = fabs( temp_jets_p4[jidx].phi() - *evt_metMuonJESCorrPhi ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_metMuonJESCorrPhi : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_metMuonJESCorrPhi );
	float dphi_tc    = fabs( temp_jets_p4[jidx].phi() - *evt_tcmetphi          ) < TMath::Pi() ? temp_jets_p4[jidx].phi() - *evt_tcmetphi          : 2*TMath::Pi() - ( temp_jets_p4[jidx].phi() - *evt_tcmetphi          );

	if( dphi < temp_dPhi_nJet_unCorrMet )
	  temp_dPhi_nJet_unCorrMet = dphi;

	if( dphi_mu < temp_dPhi_nJet_muCorrMet )
	  temp_dPhi_nJet_muCorrMet = dphi_mu;

	if( dphi_muonJESCorr < temp_dPhi_nJet_metMuonJESCorr )
	  temp_dPhi_nJet_metMuonJESCorr = dphi_muonJESCorr;

	if( dphi_tc < temp_dPhi_nJet_tcMet )
	  temp_dPhi_nJet_tcMet = dphi_tc;

	temp_sumJetPt += temp_jets_p4[jidx].pt();
      }
                                  
      hyp_dPhi_nJet_unCorrMet->push_back(temp_dPhi_nJet_unCorrMet);
      hyp_dPhi_nJet_muCorrMet->push_back(temp_dPhi_nJet_muCorrMet);
      hyp_dPhi_nJet_tcMet    ->push_back(temp_dPhi_nJet_tcMet    );   
      hyp_dPhi_nJet_metMuonJESCorr ->push_back(temp_dPhi_nJet_metMuonJESCorr );

      hyp_sumJetPt->push_back( temp_sumJetPt );

      float temp_Ht = temp_sumJetPt;

      temp_Ht += ( mus_p4->at(mus_index).pt() + els_p4->at(els_index).pt() + *evt_tcmet );

      hyp_Ht->push_back(temp_Ht);
      
      hyp_njets         ->push_back(temp_jets_p4.size()                   );     
      hyp_nojets        ->push_back(temp_other_jets_p4.size()             );     
      hyp_p4            ->push_back(mus_p4->at(mus_index)+
				    els_p4->at(els_index)                 );
	
      if(el_pt < tightptcut && mu_pt > tightptcut) {
	hyp_type            ->push_back(1);
	  
	hyp_lt_validHits    ->push_back(mus_validHits    ->at(mus_index)  );
	hyp_lt_lostHits     ->push_back(mus_lostHits     ->at(mus_index)  );
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
	hyp_lt_ptErr        ->push_back(mus_ptErr        ->at(mus_index)  );
	hyp_lt_etaErr       ->push_back(mus_etaErr       ->at(mus_index)  );
	hyp_lt_phiErr       ->push_back(mus_phiErr       ->at(mus_index)  );
	hyp_lt_p4           ->push_back(mus_p4           ->at(mus_index)  );
	hyp_lt_trk_p4       ->push_back(mus_trk_p4       ->at(mus_index)  );
	
	hyp_ll_validHits    ->push_back(els_validHits    ->at(els_index)  );
	hyp_ll_lostHits     ->push_back(els_lostHits     ->at(els_index)  );
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
	hyp_ll_ptErr        ->push_back(els_ptErr        ->at(els_index)  );
	hyp_ll_etaErr       ->push_back(els_etaErr       ->at(els_index)  );
	hyp_ll_phiErr       ->push_back(els_phiErr       ->at(els_index)  );
	hyp_ll_p4           ->push_back(els_p4           ->at(els_index)  );
	hyp_ll_trk_p4       ->push_back(els_trk_p4       ->at(els_index)  );
	
	  
      } else {
	hyp_type            ->push_back(2);
		  
	hyp_lt_validHits    ->push_back(els_validHits    ->at(els_index)  );
	hyp_lt_lostHits     ->push_back(els_lostHits     ->at(els_index)  );
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
	hyp_lt_ptErr        ->push_back(els_ptErr        ->at(els_index)  );
	hyp_lt_etaErr       ->push_back(els_etaErr       ->at(els_index)  );
	hyp_lt_phiErr       ->push_back(els_phiErr       ->at(els_index)  );
	hyp_lt_p4           ->push_back(els_p4           ->at(els_index)  );
	hyp_lt_trk_p4       ->push_back(els_trk_p4       ->at(els_index)  );
	

      
	hyp_ll_validHits    ->push_back(mus_validHits    ->at(mus_index)  );
	hyp_ll_lostHits     ->push_back(mus_lostHits     ->at(mus_index)  );
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
	hyp_ll_ptErr        ->push_back(mus_ptErr        ->at(mus_index)  );
	hyp_ll_etaErr       ->push_back(mus_etaErr       ->at(mus_index)  );
	hyp_ll_phiErr       ->push_back(mus_phiErr       ->at(mus_index)  );
	hyp_ll_p4           ->push_back(mus_p4           ->at(mus_index)  );
	hyp_ll_trk_p4       ->push_back(mus_trk_p4       ->at(mus_index)  );
	
      }
    }
  }

  iEvent.put(hyp_type                     ,"hyptype"                     );
  iEvent.put(hyp_njets                    ,"hypnjets"                    );
  iEvent.put(hyp_nojets                   ,"hypnojets"                   );
  iEvent.put(hyp_p4                      ,"hypp4"                        );
 
  iEvent.put(hyp_lt_validHits             ,"hypltvalidHits"              );
  iEvent.put(hyp_lt_lostHits              ,"hypltlostHits"               );
  iEvent.put(hyp_lt_charge                ,"hypltcharge"                 );
  iEvent.put(hyp_lt_index                 ,"hypltindex"                  );
  iEvent.put(hyp_lt_id                    ,"hypltid"                     );
  iEvent.put(hyp_lt_d0                    ,"hypltd0"                     );
  iEvent.put(hyp_lt_z0                    ,"hypltz0"                     );
  iEvent.put(hyp_lt_d0corr                ,"hypltd0corr"                 );
  iEvent.put(hyp_lt_z0corr                ,"hypltz0corr"                 );
  iEvent.put(hyp_lt_chi2                  ,"hypltchi2"                   );
  iEvent.put(hyp_lt_ndof                  ,"hypltndof"                   );
  iEvent.put(hyp_lt_d0Err                 ,"hypltd0Err"                  );
  iEvent.put(hyp_lt_z0Err                 ,"hypltz0Err"                  );
  iEvent.put(hyp_lt_ptErr                 ,"hypltptErr"                  );
  iEvent.put(hyp_lt_etaErr                ,"hypltetaErr"                 );
  iEvent.put(hyp_lt_phiErr                ,"hypltphiErr"                 );
  iEvent.put(hyp_lt_p4                    ,"hypltp4"                     );
  iEvent.put(hyp_lt_trk_p4                ,"hyplttrkp4"                  );
  
  iEvent.put(hyp_ll_validHits             ,"hypllvalidHits"              );
  iEvent.put(hyp_ll_lostHits              ,"hyplllostHits"               );
  iEvent.put(hyp_ll_charge                ,"hypllcharge"                 );
  iEvent.put(hyp_ll_index                 ,"hypllindex"                  );
  iEvent.put(hyp_ll_id                    ,"hypllid"                     );
  iEvent.put(hyp_ll_d0                    ,"hyplld0"                     );
  iEvent.put(hyp_ll_z0                    ,"hypllz0"                     );
  iEvent.put(hyp_ll_d0corr                ,"hyplld0corr"                 );
  iEvent.put(hyp_ll_z0corr                ,"hypllz0corr"                 );
  iEvent.put(hyp_ll_chi2                  ,"hypllchi2"                   );
  iEvent.put(hyp_ll_ndof                  ,"hypllndof"                   );
  iEvent.put(hyp_ll_d0Err                 ,"hyplld0Err"                  );
  iEvent.put(hyp_ll_z0Err                 ,"hypllz0Err"                  );
  iEvent.put(hyp_ll_ptErr                 ,"hypllptErr"                  );
  iEvent.put(hyp_ll_etaErr                ,"hyplletaErr"                 );
  iEvent.put(hyp_ll_phiErr                ,"hypllphiErr"                 );
  iEvent.put(hyp_ll_p4                    ,"hypllp4"                     );
  iEvent.put(hyp_ll_trk_p4                ,"hyplltrkp4"                  );
  
  iEvent.put(hyp_lt_dPhi_unCorrMet        ,"hypltdPhiunCorrMet"          );
  iEvent.put(hyp_ll_dPhi_unCorrMet        ,"hyplldPhiunCorrMet"          );
  iEvent.put(hyp_lt_dPhi_muCorrMet        ,"hypltdPhimuCorrMet"          );
  iEvent.put(hyp_ll_dPhi_muCorrMet        ,"hyplldPhimuCorrMet"          );
  iEvent.put(hyp_lt_dPhi_tcMet            ,"hypltdPhitcMet"              );
  iEvent.put(hyp_ll_dPhi_tcMet            ,"hyplldPhitcMet"              );
  iEvent.put(hyp_lt_dPhi_metMuonJESCorr    ,"hypltdPhimetMuonJESCorr"    );
  iEvent.put(hyp_ll_dPhi_metMuonJESCorr    ,"hyplldPhimetMuonJESCorr"    );
  
  iEvent.put(hyp_dPhi_nJet_unCorrMet      ,"hypdPhinJetunCorrMet"        );
  iEvent.put(hyp_dPhi_nJet_muCorrMet      ,"hypdPhinJetmuCorrMet"        );
  iEvent.put(hyp_dPhi_nJet_tcMet          ,"hypdPhinJettcMet"            );
  iEvent.put(hyp_dPhi_nJet_metMuonJESCorr ,"hypdPhinJetmetMuonJESCorr"   );
  
  iEvent.put(hyp_sumJetPt                 ,"hypsumJetPt"                 );
  iEvent.put(hyp_Ht                       ,"hypHt"                       );

  iEvent.put(hyp_mt2_tcMet                ,"hypmt2tcMet"                 );
  iEvent.put(hyp_mt2_muCorrMet            ,"hypmt2muCorrMet"             );
  iEvent.put(hyp_mt2_metMuonJESCorr       ,"hypmt2metMuonJESCorr"        );

  
  iEvent.put(hyp_jets_idx                 ,"hypjetsidx"                  );
  iEvent.put(hyp_other_jets_idx           ,"hypotherjetsidx"             );
  iEvent.put(hyp_jets_p4                  ,"hypjetsp4"                   );         
  iEvent.put(hyp_other_jets_p4            ,"hypotherjetsp4"              );       
}

// ------------ method called once each job just before starting event loop  ------------
void HypDilepMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HypDilepMaker::endJob() {
}

//----------------------------------------------------------------------------------------
bool HypDilepMaker::testJetForLeptons(const LorentzVector& jetP4, const LorentzVector& lepp4) {
  
  
  bool matched = false;
  float lepphi  = lepp4.Phi();
  float jetphi = jetP4.Phi();
   
  float lepeta  = lepp4.Eta();
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

