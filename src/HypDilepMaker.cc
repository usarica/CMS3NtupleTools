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
// $Id: HypDilepMaker.cc,v 1.9 2008/10/21 16:39:35 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/METUtilities.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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
  
  muonsInputTag            = iConfig.getParameter<InputTag>("muonsInputTag"             );
  electronsInputTag        = iConfig.getParameter<InputTag>("electronsInputTag"         );
  metInputTag              = iConfig.getParameter<InputTag>("metInputTag"               );
  jetsInputTag             = iConfig.getParameter<InputTag>("jetsInputTag"              );
  patJetsInputTag          = iConfig.getParameter<InputTag>("patJetsInputTag"            );
  trksInputTag             = iConfig.getParameter<InputTag>("trksInputTag"              );
  candToGenAssTag          = iConfig.getParameter<InputTag>("candToGenAssTag"           );  
  usingPATJets             = iConfig.getParameter<bool>    ("usingPATJets"                   );
  hypJetMinEtaCut          = iConfig.getParameter<double>  ("hypJetMinEtaCut"             );
  hypJetMaxEtaCut          = iConfig.getParameter<double>  ("hypJetMaxEtaCut"             );
  hypJetMinPtCut           = iConfig.getParameter<double>  ("hypJetMinPtCut"              );
  tightptcut               = iConfig.getParameter<double>  ("TightLepton_PtCut"         );
  looseptcut               = iConfig.getParameter<double>  ("LooseLepton_PtCut"         );



  produces<vector<int> >           ("hyptype"                  ).setBranchAlias("hyp_type"                     );
  produces<vector<int> >           ("hypnjets"                 ).setBranchAlias("hyp_njets"                    );
  produces<vector<int> >           ("hypnojets"                ).setBranchAlias("hyp_nojets"                   );  
  produces<vector<LorentzVector> > ("hypp4"                    ).setBranchAlias("hyp_p4"                       );
  
  produces<vector<int> >           ("hypltvalidHits"           ).setBranchAlias("hyp_lt_validHits"             );
  produces<vector<int> >           ("hypltlostHits"            ).setBranchAlias("hyp_lt_lostHits"              );
  produces<vector<int> >           ("hypltmcid"                ).setBranchAlias("hyp_lt_mc_id"                 );
  produces<vector<int> >           ("hypltcharge"              ).setBranchAlias("hyp_lt_charge"                );
  produces<vector<int> >           ("hypltmcmotherid"          ).setBranchAlias("hyp_lt_mc_motherid"           );
  produces<vector<int> >           ("hypltindex"               ).setBranchAlias("hyp_lt_index"                 );
  produces<vector<int> >           ("hypltid"                  ).setBranchAlias("hyp_lt_id"                    );
  produces<vector<float> >         ("hypltd0"                  ).setBranchAlias("hyp_lt_d0"                    );
  produces<vector<float> >         ("hypltz0"                  ).setBranchAlias("hyp_lt_z0"                    );
  produces<vector<float> >         ("hypltd0corr"              ).setBranchAlias("hyp_lt_d0corr"                );
  produces<vector<float> >         ("hypltz0corr"              ).setBranchAlias("hyp_lt_z0corr"                );
  produces<vector<float> >         ("hypltvertexphi"           ).setBranchAlias("hyp_lt_vertexphi"             );
  produces<vector<float> >         ("hypltchi2"                ).setBranchAlias("hyp_lt_chi2"                  );
  produces<vector<float> >         ("hypltndof"                ).setBranchAlias("hyp_lt_ndof"                  );
  produces<vector<float> >         ("hypltd0Err"               ).setBranchAlias("hyp_lt_d0Err"                 );
  produces<vector<float> >         ("hypltz0Err"               ).setBranchAlias("hyp_lt_z0Err"                 );
  produces<vector<float> >         ("hypltptErr"               ).setBranchAlias("hyp_lt_ptErr"                 );
  produces<vector<float> >         ("hypltetaErr"              ).setBranchAlias("hyp_lt_etaErr"                );
  produces<vector<float> >         ("hypltphiErr"              ).setBranchAlias("hyp_lt_phiErr"                );
  produces<vector<float> >         ("hypltouterPhi"            ).setBranchAlias("hyp_lt_outerPhi"              );
  produces<vector<float> >         ("hypltouterEta"            ).setBranchAlias("hyp_lt_outerEta"              );
  produces<vector<float> >         ("hypltiso"                 ).setBranchAlias("hyp_lt_iso"                   );
  produces<vector<float> >         ("hyplttkIso"               ).setBranchAlias("hyp_lt_tkIso"                 );
  produces<vector<LorentzVector > >("hypltp4"                  ).setBranchAlias("hyp_lt_p4"                    );
  produces<vector<LorentzVector > >("hyplttrkp4"               ).setBranchAlias("hyp_lt_trk_p4"                );
  produces<vector<LorentzVector > >("hypltmcp4"                ).setBranchAlias("hyp_lt_mc_p4"                 );
    
  
  produces<vector<int> >           ("hypllvalidHits"           ).setBranchAlias("hyp_ll_validHits"             );
  produces<vector<int> >           ("hyplllostHits"            ).setBranchAlias("hyp_ll_lostHits"              );
  produces<vector<int> >           ("hypllmcid"                ).setBranchAlias("hyp_ll_mc_id"                 );
  produces<vector<int> >           ("hypllcharge"              ).setBranchAlias("hyp_ll_charge"                );
  produces<vector<int> >           ("hypllmcmotherid"          ).setBranchAlias("hyp_ll_mc_motherid"           );
  produces<vector<int> >           ("hypllindex"               ).setBranchAlias("hyp_ll_index"                 );
  produces<vector<int> >           ("hypllid"                  ).setBranchAlias("hyp_ll_id"                    );
  produces<vector<float> >         ("hyplld0"                  ).setBranchAlias("hyp_ll_d0"                    );
  produces<vector<float> >         ("hypllz0"                  ).setBranchAlias("hyp_ll_z0"                    );
  produces<vector<float> >         ("hyplld0corr"              ).setBranchAlias("hyp_ll_d0corr"                );
  produces<vector<float> >         ("hypllz0corr"              ).setBranchAlias("hyp_ll_z0corr"                );
  produces<vector<float> >         ("hypllvertexphi"           ).setBranchAlias("hyp_ll_vertexphi"             );
  produces<vector<float> >         ("hypllchi2"                ).setBranchAlias("hyp_ll_chi2"                  );
  produces<vector<float> >         ("hypllndof"                ).setBranchAlias("hyp_ll_ndof"                  );
  produces<vector<float> >         ("hyplld0Err"               ).setBranchAlias("hyp_ll_d0Err"                 );
  produces<vector<float> >         ("hypllz0Err"               ).setBranchAlias("hyp_ll_z0Err"                 );
  produces<vector<float> >         ("hypllptErr"               ).setBranchAlias("hyp_ll_ptErr"                 );
  produces<vector<float> >         ("hyplletaErr"              ).setBranchAlias("hyp_ll_etaErr"                );
  produces<vector<float> >         ("hypllphiErr"              ).setBranchAlias("hyp_ll_phiErr"                );
  produces<vector<float> >         ("hypllouterPhi"            ).setBranchAlias("hyp_ll_outerPhi"              );
  produces<vector<float> >         ("hypllouterEta"            ).setBranchAlias("hyp_ll_outerEta"              );
  produces<vector<float> >         ("hyplliso"                 ).setBranchAlias("hyp_ll_iso"                   );
  produces<vector<float> >         ("hyplltkIso"               ).setBranchAlias("hyp_ll_tkIso"                 );
  produces<vector<LorentzVector > >("hypllp4"                  ).setBranchAlias("hyp_ll_p4"                    );
  produces<vector<LorentzVector > >("hyplltrkp4"               ).setBranchAlias("hyp_ll_trk_p4"                );
  produces<vector<LorentzVector > >("hypllmcp4"                ).setBranchAlias("hyp_ll_mc_p4"                 );



  produces<vector<float> >         ("hypmet"                   ).setBranchAlias("hyp_met"                      );
  produces<vector<float> >         ("hypmetPhi"                ).setBranchAlias("hyp_metPhi"                   );
  produces<vector<float> >         ("hypmetCaloExp"            ).setBranchAlias("hyp_metCaloExp"               );
  produces<vector<float> >         ("hypmetPhiCaloExp"         ).setBranchAlias("hyp_metPhiCaloExp"            );
  produces<vector<float> >         ("hypmetCone"               ).setBranchAlias("hyp_metCone"                  );
  produces<vector<float> >         ("hypmetPhiCone"            ).setBranchAlias("hyp_metPhiCone"               );
  produces<vector<float> >         ("hypmetNoCalo"             ).setBranchAlias("hyp_metNoCalo"                );
  produces<vector<float> >         ("hypmetPhiNoCalo"          ).setBranchAlias("hyp_metPhiNoCalo"             );
  produces<vector<float> >         ("hypmetAll"                ).setBranchAlias("hyp_metAll"                   );
  produces<vector<float> >         ("hypmetPhiAll"             ).setBranchAlias("hyp_metPhiAll"                );
  produces<vector<float> >         ("hypmetAllCaloExp"         ).setBranchAlias("hyp_metAllCaloExp"            );
  produces<vector<float> >         ("hypmetPhiAllCaloExp"      ).setBranchAlias("hyp_metPhiAllCaloExp"         );
  produces<vector<float> >         ("hypmetJes5"               ).setBranchAlias("hyp_metJes5"                  );
  produces<vector<float> >         ("hypmetPhiJes5"            ).setBranchAlias("hyp_metPhiJes5"               );
  produces<vector<float> >         ("hypmetJes10"              ).setBranchAlias("hyp_metJes10"                 );
  produces<vector<float> >         ("hypmetPhiJes10"           ).setBranchAlias("hyp_metPhiJes10"              );
  produces<vector<float> >         ("hypmetJes15"              ).setBranchAlias("hyp_metJes15"                 );
  produces<vector<float> >         ("hypmetPhiJes15"           ).setBranchAlias("hyp_metPhiJes15"              );
  produces<vector<float> >         ("hypmetJes30"              ).setBranchAlias("hyp_metJes30"                 );
  produces<vector<float> >         ("hypmetPhiJes30"           ).setBranchAlias("hyp_metPhiJes30"              );
  produces<vector<float> >         ("hypmetJes50"              ).setBranchAlias("hyp_metJes50"                 );
  produces<vector<float> >         ("hypmetPhiJes50"           ).setBranchAlias("hyp_metPhiJes50"              );
  // produces<vector<float> >         ("hypmetEMF5"               ).setBranchAlias("hyp_metEMF5"                  );
  //   produces<vector<float> >         ("hypmetPhiEMF5"            ).setBranchAlias("hyp_metPhiEMF5"               );
  //   produces<vector<float> >         ("hypmetEMF10"              ).setBranchAlias("hyp_metEMF10"                 );
  //   produces<vector<float> >         ("hypmetPhiEMF10"           ).setBranchAlias("hyp_metPhiEMF10"              );
  //   produces<vector<float> >         ("hypmetEMF15"              ).setBranchAlias("hyp_metEMF15"                 );
  //   produces<vector<float> >         ("hypmetPhiEMF15"           ).setBranchAlias("hyp_metPhiEMF15"              );
  //   produces<vector<float> >         ("hypmetEMF30"              ).setBranchAlias("hyp_metEMF30"                 );
  //   produces<vector<float> >         ("hypmetPhiEMF30"           ).setBranchAlias("hyp_metPhiEMF30"              );
  //   produces<vector<float> >         ("hypmetEMF50"              ).setBranchAlias("hyp_metEMF50"                 );
  //   produces<vector<float> >         ("hypmetPhiEMF50"           ).setBranchAlias("hyp_metPhiEMF50"              );
  produces<vector<float> >         ("hypmetDPhiJet10"          ).setBranchAlias("hyp_metDPhiJet10"             );
  produces<vector<float> >         ("hypmetDPhiJet15"          ).setBranchAlias("hyp_metDPhiJet15"             );
  produces<vector<float> >         ("hypmetDPhiJet20"          ).setBranchAlias("hyp_metDPhiJet20"             );
  produces<vector<float> >         ("hypmetDPhiTrk10"          ).setBranchAlias("hyp_metDPhiTrk10"             );
  produces<vector<float> >         ("hypmetDPhiTrk25"          ).setBranchAlias("hyp_metDPhiTrk25"             );
  produces<vector<float> >         ("hypmetDPhiTrk50"          ).setBranchAlias("hyp_metDPhiTrk50"             );

  produces<vector<vector<int> > >  ("hypjetsmcid"              ).setBranchAlias("hyp_jets_mc_id"               );
  produces<vector<vector<int> > >  ("hypotherjetsmcid"         ).setBranchAlias("hyp_other_jets_mc_id"         );
  produces<vector<vector<float> > >("hypjetsemFrac"            ).setBranchAlias("hyp_jets_emFrac"              );
  produces<vector<vector<float> > >("hypjetschFrac"            ).setBranchAlias("hyp_jets_chFrac"               );
  produces<vector<vector<float> > >("hypjetsmcemEnergy"        ).setBranchAlias("hyp_jets_mc_emEnergy"         );
  produces<vector<vector<float> > >("hypjetsmchadEnergy"       ).setBranchAlias("hyp_jets_mc_hadEnergy"        );
  produces<vector<vector<float> > >("hypjetsmcinvEnergy"       ).setBranchAlias("hyp_jets_mc_invEnergy"        );
  produces<vector<vector<float> > >("hypjetsmcotherEnergy"     ).setBranchAlias("hyp_jets_mc_otherEnergy"      );
  produces<vector<vector<float> > >("hypjetscor"               ).setBranchAlias("hyp_jets_cor"                 );
  produces<vector<vector<float> > >("hypjetsEMFcor"            ).setBranchAlias("hyp_jets_EMFcor"              );
  
  produces<vector<vector<float> > >("hypotherjetsemFrac"       ).setBranchAlias("hyp_other_jets_emFrac"        );
  produces<vector<vector<float> > >("hypotherjetschFrac"       ).setBranchAlias("hyp_other_jets_chFrac"        );
  produces<vector<vector<float> > >("hypotherjetsmcemEnergy"   ).setBranchAlias("hyp_other_jets_mc_emEnergy"   );
  produces<vector<vector<float> > >("hypotherjetsmchadEnergy"  ).setBranchAlias("hyp_other_jets_mc_hadEnergy"  );
  produces<vector<vector<float> > >("hypotherjetsmcinvEnergy"  ).setBranchAlias("hyp_other_jets_mc_invEnergy"  );
  produces<vector<vector<float> > >("hypotherjetsmcotherEnergy").setBranchAlias("hyp_other_jets_mc_otherEnergy");
  produces<vector<vector<float> > >("hypotherjetscor"          ).setBranchAlias("hyp_other_jets_cor"           );
  produces<vector<vector<float> > >("hypotherjetsEMFcor"       ).setBranchAlias("hyp_other_jets_EMFcor"        );
  
  
  produces<vector<vector<LorentzVector> > >  ("hypjetsp4"                 ).setBranchAlias("hyp_jets_p4"                     );
  produces<vector<vector<LorentzVector> > >  ("hypjetsmcp4"               ).setBranchAlias("hyp_jets_mc_p4"                  );
  produces<vector<vector<LorentzVector> > >  ("hypjetsmcgpp4"             ).setBranchAlias("hyp_jets_mc_gp_p4"               );
  produces<vector<vector<LorentzVector> > >  ("hypotherjetsp4"            ).setBranchAlias("hyp_other_jets_p4"               );
  produces<vector<vector<LorentzVector> > >  ("hypotherjetsmcp4"          ).setBranchAlias("hyp_other_jets_mc_p4"            );
  produces<vector<vector<LorentzVector> > >  ("hypotherjetsmcgpp4"        ).setBranchAlias("hyp_other_jets_mc_gp_p4"         );
  produces<vector<vector<LorentzVector> > >  ("hypotherjetspatgenPartonp4" ).setBranchAlias("hyp_other_jets_pat_genParton_p4"  );
  produces<vector<vector<LorentzVector> > >  ("hypotherjetspatgenPartonMotherp4"  ).setBranchAlias("hyp_other_jets_pat_genPartonMother_p4");
  
  if(usingPATJets) {
    produces<vector<vector<int> > >           ("hypjetspatgenPartonid"           ).setBranchAlias("hyp_jets_pat_genParton_id"             );
    produces<vector<vector<int> > >            ("hypjetspatgenPartonMotherid"     ).setBranchAlias("hyp_jets_pat_genPartonMother_id"       );
    produces<vector<vector<int> > >            ("hypjetspatpartonFlavour"         ).setBranchAlias("hyp_jets_pat_partonFlavour"            );
    produces<vector<vector<float> > >          ("hypjetspatnoCorrF"               ).setBranchAlias("hyp_jets_pat_noCorrF"                  );
    produces<vector<vector<float> > >          ("hypjetspatudsCorrF"              ).setBranchAlias("hyp_jets_pat_udsCorrF"                 );
    produces<vector<vector<float> > >          ("hypjetspatgluCorrF"              ).setBranchAlias("hyp_jets_pat_gluCorrF"                 );
    produces<vector<vector<float> > >          ("hypjetspatcCorrF"                ).setBranchAlias("hyp_jets_pat_cCorrF"                   );
    produces<vector<vector<float> > >          ("hypjetspatbCorrF"                ).setBranchAlias("hyp_jets_pat_bCorrF"                   );
    produces<vector<vector<float> > >          ("hypjetspatjetCharge"             ).setBranchAlias("hyp_jets_pat_jetCharge"                );
    produces<vector<vector<int> > >            ("hypotherjetspatgenPartonid"      ).setBranchAlias("hyp_other_jets_pat_genParton_id"       );
    produces<vector<vector<int> > >            ("hypotherjetspatgenPartonMotherid").setBranchAlias("hyp_other_jets_pat_genPartonMother_id" );
    produces<vector<vector<int> > >            ("hypotherjetspatpartonFlavour"     ).setBranchAlias("hyp_other_jets_pat_partonFlavour"      );
    produces<vector<vector<float> > >          ("hypotherjetspatnoCorrF"          ).setBranchAlias("hyp_other_jets_pat_noCorrF"            );
    produces<vector<vector<float> > >          ("hypotherjetspatudsCorrF"         ).setBranchAlias("hyp_other_jets_pat_udsCorrF"           );
    produces<vector<vector<float> > >          ("hypotherjetspatgluCorrF"         ).setBranchAlias("hyp_other_jets_pat_gluCorrF"           );
    produces<vector<vector<float> > >          ("hypotherjetspatcCorrF"           ).setBranchAlias("hyp_other_jets_pat_cCorrF"             );
    produces<vector<vector<float> > >          ("hypotherjetspatbCorrF"           ).setBranchAlias("hyp_other_jets_pat_bCorrF"             );
    produces<vector<vector<float> > >          ("hypotherjetspatjetCharge"        ).setBranchAlias("hyp_other_jets_pat_jetCharge"          );
    produces<vector<vector<LorentzVector> > >  ("hypjetspatgenPartonp4"           ).setBranchAlias("hyp_jets_pat_genParton_p4"             );
    produces<vector<vector<LorentzVector> > >  ("hypjetspatgenPartonMotherp4"     ).setBranchAlias("hyp_jets_pat_genPartonMother_p4"       );
    produces<vector<vector<LorentzVector> > >  ("hypjetspatgenJetp4"              ).setBranchAlias("hyp_jets_pat_p4"                       );
    produces<vector<vector<LorentzVector> > >  ("hypjetspatjetp4"                 ).setBranchAlias("hyp_jets_pat_jet_p4"                   );
    produces<vector<vector<LorentzVector> > >  ("hypotherjetspatgenJetp4"         ).setBranchAlias("hyp_other_jets_pat_genJet_p4"          );
    produces<vector<vector<LorentzVector> > >  ("hypotherjetspatjetp4"            ).setBranchAlias("hyp_other_jets_pat_jet_p4"             );
  }
  
  
  

    

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
  auto_ptr<vector<int> >   hyp_type                     (new vector<int>);
  auto_ptr<vector<int> >   hyp_njets                    (new vector<int>);
  auto_ptr<vector<int> >   hyp_nojets                   (new vector<int>);
  auto_ptr<vector<LorentzVector> > hyp_p4               (new vector<LorentzVector>);

  auto_ptr<vector<int> >   hyp_lt_validHits             (new vector<int>);
  auto_ptr<vector<int> >   hyp_lt_lostHits              (new vector<int>);
  auto_ptr<vector<int> >   hyp_lt_mc_id                 (new vector<int>);
  auto_ptr<vector<int> >   hyp_lt_charge                (new vector<int>);
  auto_ptr<vector<int> >   hyp_lt_mc_motherid           (new vector<int>);
  auto_ptr<vector<int> >   hyp_lt_index                 (new vector<int>);
  auto_ptr<vector<int> >   hyp_lt_id                    (new vector<int>);
  auto_ptr<vector<float> > hyp_lt_d0                    (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_z0                    (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_d0corr                (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_z0corr                (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_vertexphi             (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_chi2                  (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_ndof                  (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_d0Err                 (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_z0Err                 (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_ptErr                 (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_etaErr                (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_phiErr                (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_outerPhi              (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_outerEta              (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_iso                   (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_tkIso                 (new vector<float>);
  auto_ptr<vector<LorentzVector> > hyp_lt_p4            (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > hyp_lt_trk_p4        (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > hyp_lt_mc_p4         (new vector<LorentzVector>);

  auto_ptr<vector<int> >   hyp_ll_validHits             (new vector<int>);
  auto_ptr<vector<int> >   hyp_ll_lostHits              (new vector<int>);
  auto_ptr<vector<int> >   hyp_ll_mc_id                 (new vector<int>);
  auto_ptr<vector<int> >   hyp_ll_charge                (new vector<int>);
  auto_ptr<vector<int> >   hyp_ll_mc_motherid           (new vector<int>);
  auto_ptr<vector<int> >   hyp_ll_index                 (new vector<int>);
  auto_ptr<vector<int> >   hyp_ll_id                    (new vector<int>);
  auto_ptr<vector<float> > hyp_ll_d0                    (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_z0                    (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_d0corr                (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_z0corr                (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_vertexphi             (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_chi2                  (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_ndof                  (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_d0Err                 (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_z0Err                 (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_ptErr                 (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_etaErr                (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_phiErr                (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_outerPhi              (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_outerEta              (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_iso                   (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_tkIso                 (new vector<float>);
  auto_ptr<vector<LorentzVector> > hyp_ll_p4            (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > hyp_ll_trk_p4        (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > hyp_ll_mc_p4         (new vector<LorentzVector>);

  
  auto_ptr<vector<float> > hyp_met                      (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhi                   (new vector<float>);
  auto_ptr<vector<float> > hyp_metCaloExp               (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiCaloExp            (new vector<float>);
  auto_ptr<vector<float> > hyp_metCone                  (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiCone               (new vector<float>);
  auto_ptr<vector<float> > hyp_metNoCalo                (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiNoCalo             (new vector<float>);
  auto_ptr<vector<float> > hyp_metAll                   (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiAll                (new vector<float>);
  auto_ptr<vector<float> > hyp_metAllCaloExp            (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiAllCaloExp         (new vector<float>);
  auto_ptr<vector<float> > hyp_metJes5                  (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiJes5               (new vector<float>);
  auto_ptr<vector<float> > hyp_metJes10                 (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiJes10              (new vector<float>);
  auto_ptr<vector<float> > hyp_metJes15                 (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiJes15              (new vector<float>);
  auto_ptr<vector<float> > hyp_metJes30                 (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiJes30              (new vector<float>);
  auto_ptr<vector<float> > hyp_metJes50                 (new vector<float>);
  auto_ptr<vector<float> > hyp_metPhiJes50              (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metEMF5                  (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metPhiEMF5               (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metEMF10                 (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metPhiEMF10              (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metEMF15                 (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metPhiEMF15              (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metEMF30                 (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metPhiEMF30              (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metEMF50                 (new vector<float>);
  //   auto_ptr<vector<float> > hyp_metPhiEMF50              (new vector<float>);
  auto_ptr<vector<float> > hyp_metDPhiJet10             (new vector<float>);
  auto_ptr<vector<float> > hyp_metDPhiJet15             (new vector<float>);
  auto_ptr<vector<float> > hyp_metDPhiJet20             (new vector<float>);
  auto_ptr<vector<float> > hyp_metDPhiTrk10             (new vector<float>);
  auto_ptr<vector<float> > hyp_metDPhiTrk25             (new vector<float>);
  auto_ptr<vector<float> > hyp_metDPhiTrk50             (new vector<float>);


  auto_ptr<vector<vector<int> > >  hyp_jets_mc_id                (new vector<vector<int> >);
  auto_ptr<vector<vector<int> > >  hyp_other_jets_mc_id          (new vector<vector<int> >);

  auto_ptr<vector<vector<float> > > hyp_jets_emFrac              (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_chFrac              (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_mc_emEnergy         (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_mc_hadEnergy        (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_mc_invEnergy        (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_mc_otherEnergy      (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_cor                 (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_EMFcor              (new vector<vector<float> >);
  
  auto_ptr<vector<vector<float> > > hyp_other_jets_emFrac        (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_chFrac        (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_mc_emEnergy   (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_mc_hadEnergy  (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_mc_invEnergy  (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_mc_otherEnergy(new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_cor           (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_EMFcor        (new vector<vector<float> >);
  
  auto_ptr<vector<vector<LorentzVector> > >  hyp_jets_p4         (new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_jets_mc_p4      (new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_jets_mc_gp_p4   (new vector<vector<LorentzVector> >);

  auto_ptr<vector<vector<LorentzVector> > >  hyp_other_jets_p4      (new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_other_jets_mc_p4   (new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_other_jets_mc_gp_p4(new vector<vector<LorentzVector> >);

  auto_ptr<vector<vector<int> > >  hyp_jets_pat_genParton_id                (new vector<vector<int> >);
  auto_ptr<vector<vector<int> > >  hyp_jets_pat_genPartonMother_id          (new vector<vector<int> >);
  auto_ptr<vector<vector<int> > >  hyp_jets_pat_partonFlavour               (new vector<vector<int> >);
  auto_ptr<vector<vector<float> > > hyp_jets_pat_noCorrF                  (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_pat_udsCorrF                 (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_pat_gluCorrF                 (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_pat_cCorrF                   (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_pat_bCorrF                   (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_jets_pat_jetCharge                (new vector<vector<float> >);
  auto_ptr<vector<vector<int> > >  hyp_other_jets_pat_genParton_id          (new vector<vector<int> >);
  auto_ptr<vector<vector<int> > >  hyp_other_jets_pat_genPartonMother_id    (new vector<vector<int> >);
  auto_ptr<vector<vector<int> > >  hyp_other_jets_pat_partonFlavour         (new vector<vector<int> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_pat_noCorrF                  (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_pat_udsCorrF                 (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_pat_gluCorrF                 (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_pat_cCorrF                   (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_pat_bCorrF                   (new vector<vector<float> >);
  auto_ptr<vector<vector<float> > > hyp_other_jets_pat_jetCharge                (new vector<vector<float> >);
  
  auto_ptr<vector<vector<LorentzVector> > >  hyp_jets_pat_genParton_p4      (new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_jets_pat_genPartonMother_p4(new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_jets_pat_genJet_p4         (new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_jets_pat_jet_p4            (new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_other_jets_pat_genParton_p4(new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_other_jets_pat_genPartonMother_p4(new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_other_jets_pat_genJet_p4   (new vector<vector<LorentzVector> >);
  auto_ptr<vector<vector<LorentzVector> > >  hyp_other_jets_pat_jet_p4      (new vector<vector<LorentzVector> >);
  
  

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

  //vertex Phi
  InputTag mus_vertexphi_tag(muonsInputTag.label(),"musvertexphi");
  Handle<vector<float> > mus_vertexphi_h;
  iEvent.getByLabel(mus_vertexphi_tag, mus_vertexphi_h);
  const vector<float> *mus_vertexphi = mus_vertexphi_h.product();

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

  //outerPhi
  InputTag mus_outerPhi_tag(muonsInputTag.label(),"musouterPhi");
  Handle<vector<float> > mus_outerPhi_h;
  iEvent.getByLabel(mus_outerPhi_tag, mus_outerPhi_h);
  const vector<float> *mus_outerPhi = mus_outerPhi_h.product();


  //outerEta
  InputTag mus_outerEta_tag(muonsInputTag.label(),"musouterEta");
  Handle<vector<float> > mus_outerEta_h;
  iEvent.getByLabel(mus_outerEta_tag, mus_outerEta_h);
  const vector<float> *mus_outerEta = mus_outerEta_h.product();

  //temp iso for CMS1 validation
  InputTag mus_iso_tag(muonsInputTag.label(), "musiso");
  Handle<vector<float> > mus_iso_h;
  iEvent.getByLabel(mus_iso_tag, mus_iso_h);
  const vector<float> *mus_iso = mus_iso_h.product();

  //track isolation in dR = 0.3 from the muon object
  InputTag mus_iso03_sumPt_tag(muonsInputTag.label(),"musiso03sumPt");
  Handle<vector<float> > mus_iso03_sumPt_h;
  iEvent.getByLabel(mus_iso03_sumPt_tag, mus_iso03_sumPt_h);
  const vector<float> *mus_iso03_sumPt = mus_iso03_sumPt_h.product();

  
  //muon track P4
  InputTag mus_trk_p4_tag(muonsInputTag.label(),"mustrkp4");
  Handle<vector<LorentzVector> > mus_trk_p4_h;
  iEvent.getByLabel(mus_trk_p4_tag, mus_trk_p4_h);
  const vector<LorentzVector> *mus_trk_p4 = mus_trk_p4_h.product();

  //energy deposited in EM cal
  InputTag mus_e_em_tag(muonsInputTag.label(), "museem");
  Handle<vector<float> > mus_e_em_h;
  iEvent.getByLabel(mus_e_em_tag, mus_e_em_h);
  const vector<float> *mus_e_em = mus_e_em_h.product();

  //energy deposited in HAD cal
  InputTag mus_e_had_tag(muonsInputTag.label(), "musehad");
  Handle<vector<float> > mus_e_had_h;
  iEvent.getByLabel(mus_e_had_tag, mus_e_had_h);
  const vector<float> *mus_e_had = mus_e_had_h.product();
  
  //energy deposited in HO cal
  InputTag mus_e_ho_tag(muonsInputTag.label(), "museho");
  Handle<vector<float> > mus_e_ho_h;
  iEvent.getByLabel(mus_e_ho_tag, mus_e_ho_h);
  const vector<float> *mus_e_ho = mus_e_ho_h.product();


  //energy deposited in EM - S9
  InputTag mus_e_emS9_tag(muonsInputTag.label(), "museemS9");
  Handle<vector<float> > mus_e_emS9_h;
  iEvent.getByLabel(mus_e_emS9_tag, mus_e_emS9_h);
  const vector<float> *mus_e_emS9 = mus_e_emS9_h.product();

  //energy deposited in HAD - s9
  InputTag mus_e_hadS9_tag(muonsInputTag.label(), "musehadS9");
  Handle<vector<float> > mus_e_hadS9_h;
  iEvent.getByLabel(mus_e_hadS9_tag, mus_e_hadS9_h);
  const vector<float> *mus_e_hadS9 = mus_e_hadS9_h.product();
  
  //energy deposited in HO cal
  InputTag mus_e_hoS9_tag(muonsInputTag.label(), "musehoS9");
  Handle<vector<float> > mus_e_hoS9_h;
  iEvent.getByLabel(mus_e_hoS9_tag, mus_e_hoS9_h);
  const vector<float> *mus_e_hoS9 = mus_e_hoS9_h.product();  



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

  
  //vertex Phi
  InputTag els_vertexphi_tag(electronsInputTag.label(),"elsvertexphi");
  Handle<vector<float> > els_vertexphi_h;
  iEvent.getByLabel(els_vertexphi_tag, els_vertexphi_h);
  const vector<float> *els_vertexphi = els_vertexphi_h.product();

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

  //outerPhi
  InputTag els_outerPhi_tag(electronsInputTag.label(),"elsouterPhi");
  Handle<vector<float> > els_outerPhi_h;
  iEvent.getByLabel(els_outerPhi_tag, els_outerPhi_h);
  const vector<float> *els_outerPhi = els_outerPhi_h.product();


  //outerEta
  InputTag els_outerEta_tag(electronsInputTag.label(),"elsouterEta");
  Handle<vector<float> > els_outerEta_h;
  iEvent.getByLabel(els_outerEta_tag, els_outerEta_h);
  const vector<float> *els_outerEta = els_outerEta_h.product();

  //track isolation in dR = 0.3
  InputTag els_tkIso_tag(electronsInputTag.label(),"elstkIso");
  Handle<vector<float> > els_tkIso_h;
  iEvent.getByLabel(els_tkIso_tag, els_tkIso_h);
  const vector<float> *els_tkIso = els_tkIso_h.product();

  //electron track P4
  InputTag els_trk_p4_tag(electronsInputTag.label(),"elstrkp4");
  Handle<vector<LorentzVector> > els_trk_p4_h;
  iEvent.getByLabel(els_trk_p4_tag, els_trk_p4_h);
  const vector<LorentzVector> *els_trk_p4 = els_trk_p4_h.product();

  

  //--------------------------------------------------------------------
  //Get the Jet collections
  //--------------------------------------------------------------------
    
  //Get the jet EM fraction
  InputTag jets_emFrac_tag(jetsInputTag.label(), "jetsemFrac");
  Handle<vector<float> > jets_emFrac_h;
  iEvent.getByLabel(jets_emFrac_tag, jets_emFrac_h);
  const vector<float> *jets_emFrac = jets_emFrac_h.product();
  
  //Get the jet charge fraction
  InputTag jets_chFrac_tag(jetsInputTag.label(), "jetschFrac");
  Handle<vector<float> > jets_chFrac_h;
  iEvent.getByLabel(jets_chFrac_tag, jets_chFrac_h);
  const vector<float> *jets_chFrac = jets_chFrac_h.product();

  //Get the JES correction collection
  InputTag jetscor_tag(jetsInputTag.label(), "jetscor");
  Handle<vector<float> > jetscor_h;
  iEvent.getByLabel(jetscor_tag, jetscor_h);
  const vector<float> *jets_cor = jetscor_h.product();

  //Get the EMF correction collection
  InputTag jetsEMF_tag(jetsInputTag.label(), "jetsEMFcor");
  Handle<vector<float> > jetsEMF_h;
  iEvent.getByLabel(jetsEMF_tag, jetsEMF_h);
  const vector<float> *jets_EMFcor = jetsEMF_h.product();
  
  //jet p4
  InputTag jets_p4_tag(jetsInputTag.label(), "jetsp4");
  Handle<vector<LorentzVector> > jets_p4_h;
  iEvent.getByLabel(jets_p4_tag, jets_p4_h);
  const vector<LorentzVector> *jets_p4 = jets_p4_h.product();

  //------------------------------------------------------------
  //Get the PAT jet corrections if we're using PATAF Jets
  //-----------------------------------------------------------
  InputTag jets_pat_genParton_id_tag(patJetsInputTag.label(), "jetspatgenPartonid");
  Handle<vector<int> > jets_pat_genParton_id_h;
    
  InputTag jets_pat_genPartonMother_id_tag(patJetsInputTag.label(), "jetspatgenPartonMotherid");
  Handle<vector<int> > jets_pat_genPartonMother_id_h;
    
  InputTag jets_pat_partonFlavour_tag(patJetsInputTag.label(), "jetspatpartonFlavour");
  Handle<vector<int> > jets_pat_partonFlavour_h;

  
  
  InputTag jets_pat_genParton_p4_tag(patJetsInputTag.label(), "jetspatgenPartonp4");
  Handle<vector<LorentzVector> > jets_pat_genParton_p4_h;
    
  InputTag jets_pat_genPartonMother_p4_tag(patJetsInputTag.label(), "jetspatgenPartonMotherp4");
  Handle<vector<LorentzVector> > jets_pat_genPartonMother_p4_h;
    
  InputTag jets_pat_genJet_p4_tag(patJetsInputTag.label(), "jetspatgenJetp4");
  Handle<vector<LorentzVector> > jets_pat_genJet_p4_h;
    
  //This is corrected! Be careful when using this to correct the MET!
  InputTag jets_pat_jet_p4_tag(patJetsInputTag.label(), "jetspatjetp4");
  Handle<vector<LorentzVector> > jets_pat_jet_p4_h;
    
  //correction factor
  InputTag jets_pat_noCorrF_tag(patJetsInputTag.label(), "jetspatnoCorrF");
  Handle<vector<float> > jets_pat_noCorrF_h;
  
  InputTag jets_pat_udsCorrF_tag(patJetsInputTag.label(),"jetspatudsCorrF");
  Handle<vector<float> > jets_pat_udsCorrF_h;

  InputTag jets_pat_gluCorrF_tag(patJetsInputTag.label(),"jetspatgluCorrF");
  Handle<vector<float> > jets_pat_gluCorrF_h;

  InputTag jets_pat_cCorrF_tag(patJetsInputTag.label(),"jetspatcCorrF");
  Handle<vector<float> > jets_pat_cCorrF_h;

  InputTag jets_pat_bCorrF_tag(patJetsInputTag.label(),"jetspatbCorrF");
  Handle<vector<float> > jets_pat_bCorrF_h;

  InputTag jets_pat_jetCharge_tag(patJetsInputTag.label(),"jetspatjetCharge");
  Handle<vector<float> > jets_pat_jetCharge_h;
  

    
  if(usingPATJets) {
    iEvent.getByLabel(jets_pat_genParton_id_tag, jets_pat_genParton_id_h);
    
    iEvent.getByLabel(jets_pat_genPartonMother_id_tag, jets_pat_genPartonMother_id_h);
    
    iEvent.getByLabel(jets_pat_partonFlavour_tag, jets_pat_partonFlavour_h);
    
    iEvent.getByLabel(jets_pat_genParton_p4_tag, jets_pat_genParton_p4_h);
            
    iEvent.getByLabel(jets_pat_genPartonMother_p4_tag, jets_pat_genPartonMother_p4_h);
        
    iEvent.getByLabel(jets_pat_genJet_p4_tag, jets_pat_genJet_p4_h);
        
    iEvent.getByLabel(jets_pat_jet_p4_tag, jets_pat_jet_p4_h);
        
    iEvent.getByLabel(jets_pat_noCorrF_tag, jets_pat_noCorrF_h);

    iEvent.getByLabel(jets_pat_udsCorrF_tag, jets_pat_udsCorrF_h);
        
    iEvent.getByLabel(jets_pat_gluCorrF_tag,jets_pat_gluCorrF_h);
        
    iEvent.getByLabel(jets_pat_cCorrF_tag, jets_pat_cCorrF_h);
        
    iEvent.getByLabel(jets_pat_bCorrF_tag, jets_pat_bCorrF_h);
        
    iEvent.getByLabel(jets_pat_jetCharge_tag,jets_pat_jetCharge_h);
            
  }
 
  //event met - this is uncorrected
  InputTag met_tag(metInputTag.label(), "evtmet");
  Handle<float> met_tag_h;
  iEvent.getByLabel(met_tag, met_tag_h);
  const float* evt_met = met_tag_h.product();

  //event metPhi
  InputTag metphi_tag(metInputTag.label(), "evtmetPhi");
  Handle<float> metphi_tag_h;
  iEvent.getByLabel(metphi_tag, metphi_tag_h);
  const float* evt_metphi = metphi_tag_h.product();

  //track p4s
  InputTag trks_p4_tag(trksInputTag.label(), "trkstrkp4");
  Handle<vector<LorentzVector> > trks_p4_h;
  iEvent.getByLabel(trks_p4_tag, trks_p4_h);
  const vector<LorentzVector>* trks_p4 = trks_p4_h.product();



  //Generator to Candidate matching stuff

  //PDG id of matched MC particle
  InputTag mus_mc_id_tag(candToGenAssTag.label(),"musmcid");
  Handle<vector<int> > mus_mc_id_h;
  iEvent.getByLabel(mus_mc_id_tag, mus_mc_id_h);
  const vector<int> *mus_mc_id = mus_mc_id_h.product();

   //PDG id of MC matched mother 
  InputTag mus_mc_motherid_tag(candToGenAssTag.label(),"musmcmotherid");
  Handle<vector<int> > mus_mc_motherid_h;
  iEvent.getByLabel(mus_mc_motherid_tag, mus_mc_motherid_h);
  const vector<int> *mus_mc_motherid = mus_mc_motherid_h.product();

  //muon mc P4
  InputTag mus_mc_p4_tag(candToGenAssTag.label(),"musmcp4");
  Handle<vector<LorentzVector> > mus_mc_p4_h;
  iEvent.getByLabel(mus_mc_p4_tag, mus_mc_p4_h);
  const vector<LorentzVector> *mus_mc_p4 = mus_mc_p4_h.product();

  //PDG id of matched MC particle
  InputTag els_mc_id_tag(candToGenAssTag.label(),"elsmcid");
  Handle<vector<int> > els_mc_id_h;
  iEvent.getByLabel(els_mc_id_tag, els_mc_id_h);
  const vector<int> *els_mc_id = els_mc_id_h.product();

  //PDG id of MC matched mother 
  InputTag els_mc_motherid_tag(candToGenAssTag.label(),"elsmcmotherid");
  Handle<vector<int> > els_mc_motherid_h;
  iEvent.getByLabel(els_mc_motherid_tag, els_mc_motherid_h);
  const vector<int> *els_mc_motherid = els_mc_motherid_h.product();

  
  //electron mc P4
  InputTag els_mc_p4_tag(candToGenAssTag.label(),"elsmcp4");
  Handle<vector<LorentzVector> > els_mc_p4_h;
  iEvent.getByLabel(els_mc_p4_tag, els_mc_p4_h);
  const vector<LorentzVector> *els_mc_p4 = els_mc_p4_h.product();

  //jet MC matching
  InputTag jets_mc_id_tag(candToGenAssTag.label(), "jetsmcid");
  Handle<vector<int> > jets_mc_id_h;
  iEvent.getByLabel(jets_mc_id_tag, jets_mc_id_h);
  const vector<int> *jets_mc_id = jets_mc_id_h.product();

  //Get the jet mc EM energy
  InputTag jets_mc_emEnergy_tag(candToGenAssTag.label(), "jetsmcemEnergy");
  Handle<vector<float> > jets_mc_emEnergy_h;
  iEvent.getByLabel(jets_mc_emEnergy_tag, jets_mc_emEnergy_h);
  const vector<float> *jets_mc_emEnergy = jets_mc_emEnergy_h.product();

  //Get the jet mc had energy
  InputTag jets_mc_hadEnergy_tag(candToGenAssTag.label(), "jetsmchadEnergy");
  Handle<vector<float> > jets_mc_hadEnergy_h;
  iEvent.getByLabel(jets_mc_hadEnergy_tag, jets_mc_hadEnergy_h);
  const vector<float> *jets_mc_hadEnergy = jets_mc_hadEnergy_h.product();

  //Get the jet mc invisible energy
  InputTag jets_mc_invEnergy_tag(candToGenAssTag.label(), "jetsmcinvEnergy");
  Handle<vector<float> > jets_mc_invEnergy_h;
  iEvent.getByLabel(jets_mc_invEnergy_tag, jets_mc_invEnergy_h);
  const vector<float> *jets_mc_invEnergy = jets_mc_invEnergy_h.product();

  //Get the jet mc other energy
  InputTag jets_mc_otherEnergy_tag(candToGenAssTag.label(), "jetsmcotherEnergy");
  Handle<vector<float> > jets_mc_otherEnergy_h;
  iEvent.getByLabel(jets_mc_otherEnergy_tag, jets_mc_otherEnergy_h);
  const vector<float> *jets_mc_otherEnergy = jets_mc_otherEnergy_h.product();

  //jet mc p4
  InputTag jets_mc_p4_tag(candToGenAssTag.label(), "jetsmcp4");
  Handle<vector<LorentzVector> > jets_mc_p4_h;
  iEvent.getByLabel(jets_mc_p4_tag, jets_mc_p4_h);
  const vector<LorentzVector> *jets_mc_p4 = jets_mc_p4_h.product();

  //jet mc genParton p4
  InputTag jets_mc_gp_p4_tag(candToGenAssTag.label(), "jetsmcgpp4");
  Handle<vector<LorentzVector> > jets_mc_gp_p4_h;
  iEvent.getByLabel(jets_mc_gp_p4_tag, jets_mc_gp_p4_h);
  const vector<LorentzVector> *jets_mc_gp_p4 = jets_mc_gp_p4_h.product();
  


  double metAll           = *evt_met;
  double metPhiAll        = *evt_metphi;
  double metAllCaloExp    = *evt_met;
  double metPhiAllCaloExp = *evt_metphi;

  
  //calculate the MET corrections that use all muons here
  //since its the same for all the candidate pairs
  for(unsigned int i = 0; i < mus_p4->size(); i++) {

    pair<LorentzVector, LorentzVector> muon_pair = make_pair(mus_p4->at(i),
							     mus_trk_p4->at(i) );
    METUtilities::correctMETmuons_crossedE(muon_pair, metAll, metPhiAll,
					   mus_e_em->at(i), mus_e_had->at(i),  mus_e_ho->at(i) );
    METUtilities::correctMETmuons_expMIP(muon_pair, metAllCaloExp,
					 metPhiAllCaloExp );
  }

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
      vector<int>            temp_jets_mc_id;
      vector<float>          temp_jets_emFrac;
      vector<float>          temp_jets_chFrac;
      vector<float>          temp_jets_mc_emEnergy;
      vector<float>          temp_jets_mc_hadEnergy;
      vector<float>          temp_jets_mc_invEnergy;
      vector<float>          temp_jets_mc_otherEnergy;
      vector<float>          temp_jets_cor;
      vector<float>          temp_jets_EMFcor;
                                                   
      vector<LorentzVector>  temp_jets_p4;
      vector<LorentzVector>  temp_jets_mc_p4;
      vector<LorentzVector>  temp_jets_mc_gp_p4;

      vector<int>            temp_other_jets_mc_id;
      vector<float>          temp_other_jets_emFrac;
      vector<float>          temp_other_jets_chFrac;
      vector<float>          temp_other_jets_mc_emEnergy;
      vector<float>          temp_other_jets_mc_hadEnergy;
      vector<float>          temp_other_jets_mc_invEnergy;
      vector<float>          temp_other_jets_mc_otherEnergy;
      vector<float>          temp_other_jets_cor;
      vector<float>          temp_other_jets_EMFcor;
      
      vector<LorentzVector>  temp_other_jets_p4;
      vector<LorentzVector>  temp_other_jets_mc_p4;
      vector<LorentzVector>  temp_other_jets_mc_gp_p4;

      vector<int>            temp_jets_pat_genParton_id;
      vector<int>            temp_jets_pat_genPartonMother_id;
      vector<int>            temp_jets_pat_partonFlavour;
      vector<float>          temp_jets_pat_noCorrF;                  
      vector<float>          temp_jets_pat_udsCorrF;                 
      vector<float>          temp_jets_pat_gluCorrF;                 
      vector<float>          temp_jets_pat_cCorrF;                   
      vector<float>          temp_jets_pat_bCorrF;                   
      vector<float>          temp_jets_pat_jetCharge;                
      vector<LorentzVector>  temp_jets_pat_genParton_p4;
      vector<LorentzVector>  temp_jets_pat_genPartonMother_p4;
      vector<LorentzVector>  temp_jets_pat_genJet_p4;
      vector<LorentzVector>  temp_jets_pat_jet_p4;
      
      vector<int>            temp_other_jets_pat_genParton_id;
      vector<int>            temp_other_jets_pat_genPartonMother_id;
      vector<int>            temp_other_jets_pat_partonFlavour;
      vector<float>          temp_other_jets_pat_noCorrF;                  
      vector<float>          temp_other_jets_pat_udsCorrF;                 
      vector<float>          temp_other_jets_pat_gluCorrF;                 
      vector<float>          temp_other_jets_pat_cCorrF;                   
      vector<float>          temp_other_jets_pat_bCorrF;                   
      vector<float>          temp_other_jets_pat_jetCharge;                
      vector<LorentzVector>  temp_other_jets_pat_genParton_p4;
      vector<LorentzVector>  temp_other_jets_pat_genPartonMother_p4;
      vector<LorentzVector>  temp_other_jets_pat_genJet_p4;
      vector<LorentzVector>  temp_other_jets_pat_jet_p4;
      	
      
      //these are for the MET correction later
      vector<LorentzVector> jets_noel_p4;
      vector<float> jets_noel_jescor;
      //vector<float> jets_noel_EMFcor;
	
      for(unsigned int i = 0; i<jets_p4->size(); i++) {

	//These have to uncorrected, so if we use PATjets
	//we have to scale accordingly
	if(usingPATJets) {
	  float px = (jets_p4->at(i).px())*jets_pat_noCorrF_h->at(i);
	  float py = (jets_p4->at(i).py())*jets_pat_noCorrF_h->at(i);
	  float pz = (jets_p4->at(i).pz())*jets_pat_noCorrF_h->at(i);
	  float E  = (jets_p4->at(i).E())*jets_pat_noCorrF_h->at(i);
	  LorentzVector temp(px, py, pz, E);
	  jets_noel_p4        .push_back(temp);
	  jets_noel_jescor    .push_back(jets_cor->at(i)/(jets_pat_noCorrF_h->at(i)));
	} else {
	  jets_noel_p4        .push_back(jets_p4    ->at(i));
	  jets_noel_jescor    .push_back(jets_cor   ->at(i));
	}
	
	//jets_noel_EMFcor    .push_back(jets_EMFcor->at(i));
	
	double jet_eta = jets_p4->at(i).eta();
	double jet_pt  = jets_p4->at(i).Pt();
	double jet_ptcut;
	
	if(usingPATJets) {
	  jet_ptcut =  hypJetMinPtCut/(jets_pat_noCorrF_h->at(i));
	} else jet_ptcut = hypJetMinPtCut;
	
	if( hypJetMinEtaCut < jet_eta && 
	    jet_eta < hypJetMaxEtaCut && 
	    jet_pt  > jet_ptcut) { //hyp jets
	  
	  temp_jets_mc_id                  .push_back(jets_mc_id           ->at(i));
	  temp_jets_emFrac                 .push_back(jets_emFrac          ->at(i));
	  temp_jets_chFrac                 .push_back(jets_chFrac          ->at(i));
	  
	  temp_jets_mc_emEnergy            .push_back(jets_mc_emEnergy     ->at(i));
	  temp_jets_mc_hadEnergy           .push_back(jets_mc_hadEnergy    ->at(i));
	  temp_jets_mc_invEnergy           .push_back(jets_mc_invEnergy    ->at(i));
	  temp_jets_mc_otherEnergy         .push_back(jets_mc_otherEnergy  ->at(i));  
	  temp_jets_cor                    .push_back(jets_cor             ->at(i));
	  temp_jets_EMFcor                 .push_back(jets_EMFcor          ->at(i));
	  
	  temp_jets_p4                     .push_back(jets_p4              ->at(i));
	  temp_jets_mc_p4                  .push_back(jets_mc_p4           ->at(i));
	  temp_jets_mc_gp_p4               .push_back(jets_mc_gp_p4        ->at(i));

	  
	  
	  if(usingPATJets) {
	    temp_jets_pat_genParton_id      .push_back(jets_pat_genParton_id_h       ->at(i));
	    temp_jets_pat_genPartonMother_id.push_back(jets_pat_genPartonMother_id_h ->at(i));
	    temp_jets_pat_partonFlavour     .push_back(jets_pat_partonFlavour_h      ->at(i));
	    temp_jets_pat_genParton_p4      .push_back(jets_pat_genParton_p4_h       ->at(i));
	    temp_jets_pat_genPartonMother_p4.push_back(jets_pat_genPartonMother_p4_h ->at(i));
	    temp_jets_pat_genJet_p4         .push_back(jets_pat_genJet_p4_h          ->at(i));
	    temp_jets_pat_jet_p4            .push_back(jets_pat_jet_p4_h             ->at(i));
	    temp_jets_pat_noCorrF           .push_back(jets_pat_noCorrF_h            ->at(i));                  
	    temp_jets_pat_udsCorrF          .push_back(jets_pat_udsCorrF_h           ->at(i));                 
	    temp_jets_pat_gluCorrF          .push_back(jets_pat_gluCorrF_h           ->at(i));                 
	    temp_jets_pat_cCorrF            .push_back(jets_pat_cCorrF_h             ->at(i));                   
	    temp_jets_pat_bCorrF            .push_back(jets_pat_bCorrF_h             ->at(i));                   
	    temp_jets_pat_jetCharge         .push_back(jets_pat_jetCharge_h          ->at(i));               
	  }
	  
	} else {
	  temp_other_jets_mc_id            .push_back(jets_mc_id           ->at(i));
	  temp_other_jets_emFrac           .push_back(jets_emFrac          ->at(i));
	  temp_other_jets_chFrac           .push_back(jets_chFrac          ->at(i));
	  
	  temp_other_jets_mc_emEnergy      .push_back(jets_mc_emEnergy     ->at(i));
	  temp_other_jets_mc_hadEnergy     .push_back(jets_mc_hadEnergy    ->at(i));
	  temp_other_jets_mc_invEnergy     .push_back(jets_mc_invEnergy    ->at(i));
	  temp_other_jets_mc_otherEnergy   .push_back(jets_mc_otherEnergy  ->at(i));  
	  temp_other_jets_cor              .push_back(jets_cor             ->at(i));
	  temp_other_jets_EMFcor           .push_back(jets_EMFcor          ->at(i));
	  
	  temp_other_jets_p4               .push_back(jets_p4              ->at(i));
	  temp_other_jets_mc_p4            .push_back(jets_mc_p4           ->at(i));
	  temp_other_jets_mc_gp_p4         .push_back(jets_mc_gp_p4        ->at(i));

	   
	  if(usingPATJets) {
	    temp_other_jets_pat_genParton_id      .push_back(jets_pat_genParton_id_h       ->at(i));
	    temp_other_jets_pat_genPartonMother_id.push_back(jets_pat_genPartonMother_id_h ->at(i));
	    temp_other_jets_pat_partonFlavour     .push_back(jets_pat_partonFlavour_h      ->at(i));
	    temp_other_jets_pat_genParton_p4      .push_back(jets_pat_genParton_p4_h       ->at(i));
	    temp_other_jets_pat_genPartonMother_p4.push_back(jets_pat_genPartonMother_p4_h ->at(i));
	    temp_other_jets_pat_genJet_p4         .push_back(jets_pat_genJet_p4_h          ->at(i));
	    temp_other_jets_pat_jet_p4            .push_back(jets_pat_jet_p4_h             ->at(i));
	    temp_other_jets_pat_noCorrF           .push_back(jets_pat_noCorrF_h            ->at(i));                  
	    temp_other_jets_pat_udsCorrF          .push_back(jets_pat_udsCorrF_h           ->at(i));                 
	    temp_other_jets_pat_gluCorrF          .push_back(jets_pat_gluCorrF_h           ->at(i));                 
	    temp_other_jets_pat_cCorrF            .push_back(jets_pat_cCorrF_h             ->at(i));                   
	    temp_other_jets_pat_bCorrF            .push_back(jets_pat_bCorrF_h             ->at(i));                   
	    temp_other_jets_pat_jetCharge        .push_back(jets_pat_jetCharge_h          ->at(i));               
	  }

	}
      }
      
      //push these into the hyp_jets and hyp_other_jets vars
      hyp_jets_mc_id                 ->push_back(temp_jets_mc_id              );
      hyp_jets_emFrac                ->push_back(temp_jets_emFrac             );
      hyp_jets_chFrac                ->push_back(temp_jets_chFrac             );
      hyp_jets_mc_emEnergy           ->push_back(temp_jets_mc_emEnergy        );
      hyp_jets_mc_hadEnergy          ->push_back(temp_jets_mc_hadEnergy       );
      hyp_jets_mc_invEnergy          ->push_back(temp_jets_mc_invEnergy       );
      hyp_jets_mc_otherEnergy        ->push_back(temp_jets_mc_otherEnergy     );
      hyp_jets_cor                   ->push_back(temp_jets_cor                );
      hyp_jets_EMFcor                ->push_back(temp_jets_EMFcor             );
      
      hyp_jets_p4                    ->push_back(temp_jets_p4                  );
      hyp_jets_mc_p4                 ->push_back(temp_jets_mc_p4               );
      hyp_jets_mc_gp_p4              ->push_back(temp_jets_mc_gp_p4            );
      
      hyp_other_jets_mc_id           ->push_back(temp_other_jets_mc_id         );
      hyp_other_jets_emFrac          ->push_back(temp_other_jets_emFrac        );
      hyp_other_jets_chFrac          ->push_back(temp_other_jets_chFrac        );
      hyp_other_jets_mc_emEnergy     ->push_back(temp_other_jets_mc_emEnergy   );
      hyp_other_jets_mc_hadEnergy    ->push_back(temp_other_jets_mc_hadEnergy  );
      hyp_other_jets_mc_invEnergy    ->push_back(temp_other_jets_mc_invEnergy  );
      hyp_other_jets_mc_otherEnergy  ->push_back(temp_other_jets_mc_otherEnergy);
      hyp_other_jets_cor             ->push_back(temp_other_jets_cor           );
      hyp_other_jets_EMFcor          ->push_back(temp_other_jets_EMFcor        );
      
      hyp_other_jets_p4              ->push_back(temp_other_jets_p4            );
      hyp_other_jets_mc_p4           ->push_back(temp_other_jets_mc_p4         );
      hyp_other_jets_mc_gp_p4        ->push_back(temp_other_jets_mc_gp_p4      );
      
      
      if(usingPATJets) {
	hyp_jets_pat_genParton_id            ->push_back(temp_jets_pat_genParton_id             );
	hyp_jets_pat_genPartonMother_id      ->push_back(temp_jets_pat_genPartonMother_id       );
	hyp_jets_pat_partonFlavour           ->push_back(temp_jets_pat_partonFlavour            );
	hyp_jets_pat_genParton_p4            ->push_back(temp_jets_pat_genParton_p4             );
	hyp_jets_pat_genPartonMother_p4      ->push_back(temp_jets_pat_genPartonMother_p4       );
	hyp_jets_pat_genJet_p4               ->push_back(temp_jets_pat_genJet_p4                );
	hyp_jets_pat_jet_p4                  ->push_back(temp_jets_pat_jet_p4                   );
	hyp_jets_pat_noCorrF                 ->push_back(temp_jets_pat_noCorrF                  );                  
	hyp_jets_pat_udsCorrF                ->push_back(temp_jets_pat_udsCorrF                 );                 
	hyp_jets_pat_gluCorrF                ->push_back(temp_jets_pat_gluCorrF                 );                 
	hyp_jets_pat_cCorrF                  ->push_back(temp_jets_pat_cCorrF                   );                   
	hyp_jets_pat_bCorrF                  ->push_back(temp_jets_pat_bCorrF                   );                   
	hyp_jets_pat_jetCharge              ->push_back(temp_jets_pat_jetCharge                );               

	hyp_other_jets_pat_genParton_id      ->push_back(temp_other_jets_pat_genParton_id       );
	hyp_other_jets_pat_genPartonMother_id->push_back(temp_other_jets_pat_genPartonMother_id );
	hyp_other_jets_pat_partonFlavour     ->push_back(temp_other_jets_pat_partonFlavour      );
	hyp_other_jets_pat_genParton_p4      ->push_back(temp_other_jets_pat_genParton_p4       );
	hyp_other_jets_pat_genPartonMother_p4->push_back(temp_other_jets_pat_genPartonMother_p4 );
	hyp_other_jets_pat_genJet_p4         ->push_back(temp_other_jets_pat_genJet_p4          );
	hyp_other_jets_pat_jet_p4            ->push_back(temp_other_jets_pat_jet_p4             );
	hyp_other_jets_pat_noCorrF           ->push_back(temp_other_jets_pat_noCorrF            );                  
	hyp_other_jets_pat_udsCorrF          ->push_back(temp_other_jets_pat_udsCorrF           );                 
	hyp_other_jets_pat_gluCorrF          ->push_back(temp_other_jets_pat_gluCorrF           );                 
	hyp_other_jets_pat_cCorrF            ->push_back(temp_other_jets_pat_cCorrF             );                   
	hyp_other_jets_pat_bCorrF            ->push_back(temp_other_jets_pat_bCorrF             );                   
	hyp_other_jets_pat_jetCharge        ->push_back(temp_other_jets_pat_jetCharge          );               
      }

      //correct the met for the hyp muons
      pair<LorentzVector, LorentzVector> tightmuon_pair = make_pair(mus_p4->at(tight_index),
								    mus_trk_p4->at(tight_index) );
      pair<LorentzVector, LorentzVector> loosemuon_pair = make_pair(mus_p4->at(loose_index),
								    mus_trk_p4->at(loose_index) );
      
      double hypmet = *evt_met;
      double hypmetPhi = *evt_metphi;
      METUtilities::correctMETmuons_crossedE(tightmuon_pair,
					     hypmet, hypmetPhi, mus_e_em->at(tight_index), 
					     mus_e_had->at(tight_index),  mus_e_ho->at(tight_index) );
      METUtilities::correctMETmuons_crossedE(loosemuon_pair,
					     hypmet, hypmetPhi, mus_e_em->at(loose_index), 
					     mus_e_had->at(loose_index), mus_e_ho->at(loose_index) );
      //now use the expected MIP deposit
      double hypmet_MIP    = *evt_met;
      double hypmetPhi_MIP = *evt_metphi;
      METUtilities::correctMETmuons_expMIP(tightmuon_pair,
					   hypmet_MIP, hypmetPhi_MIP );
      METUtilities::correctMETmuons_expMIP(loosemuon_pair,
					   hypmet_MIP, hypmetPhi_MIP );
      //now the s9 correction
      double hypmet_S9    = *evt_met;
      double hypmetPhi_S9 = *evt_metphi;
      METUtilities::correctMETmuons_S9E(tightmuon_pair,
					hypmet_S9, hypmetPhi_S9, mus_e_emS9->at(tight_index), 
					mus_e_hadS9->at(tight_index),  mus_e_hoS9->at(tight_index) );
      METUtilities::correctMETmuons_S9E(loosemuon_pair,
					hypmet_S9, hypmetPhi_S9, mus_e_emS9->at(loose_index), 
					mus_e_hadS9->at(loose_index), mus_e_hoS9->at(loose_index) );
      //no calo or MIP correction
      double hypmet_nocalo = *evt_met;
      double hypmetPhi_nocalo = *evt_metphi;
      METUtilities::correctMETmuons_nocalo(tightmuon_pair,
					   hypmet_nocalo, hypmetPhi_nocalo);
      METUtilities::correctMETmuons_nocalo(loosemuon_pair,
					   hypmet_nocalo, hypmetPhi_nocalo);
      
      
      double metJes5 = hypmet;
      double metPhiJes5 = hypmetPhi;
      double metJes10 = hypmet;
      double metPhiJes10 = hypmetPhi;
      double metJes15 = hypmet;
      double metPhiJes15 = hypmetPhi;
      double metJes30 = hypmet;
      double metPhiJes30 = hypmetPhi;
      double metJes50 = hypmet;
      double metPhiJes50 = hypmetPhi;
      
      // double metEMF5 = hypmet;
      //       double metPhiEMF5 = hypmetPhi;
      //       double metEMF10 = hypmet;
      //       double metPhiEMF10 = hypmetPhi;
      //       double metEMF15 = hypmet;
      //       double metPhiEMF15 = hypmetPhi;
      //       double metEMF30 = hypmet;
      //       double metPhiEMF30 = hypmetPhi;
      //       double metEMF50 = hypmet;
      //       double metPhiEMF50 = hypmetPhi;
      
      
      METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				    metJes5, metPhiJes5, 5);
      METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				    metJes10, metPhiJes10, 10);
      METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				    metJes15, metPhiJes15, 15);
      METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				    metJes30, metPhiJes30, 30);
      METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				    metJes50, metPhiJes50, 50);
      
             
      //   METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
      // 				    metEMF5, metPhiEMF5, 5);
      //       METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
      // 				    metEMF10, metPhiEMF10, 10);
      //       METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
      // 				    metEMF15, metPhiEMF15, 15);
      //       METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
      // 				    metEMF30, metPhiEMF30, 30);
      //       METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
      // 				    metEMF50, metPhiEMF50, 50);
      
      

      hyp_type          ->push_back(0                                     );
      hyp_njets         ->push_back(temp_jets_p4.size()                   );     
      hyp_nojets        ->push_back(temp_other_jets_p4.size()             );     
      hyp_p4            ->push_back(mus_p4->at(tight_index)+
				    mus_p4->at(loose_index)               );
      
      hyp_lt_validHits    ->push_back(mus_validHits    ->at(tight_index)  );
      hyp_lt_lostHits     ->push_back(mus_lostHits     ->at(tight_index)  );
      hyp_lt_mc_id        ->push_back(mus_mc_id        ->at(tight_index)  );
      hyp_lt_charge       ->push_back(mus_charge       ->at(tight_index)  );
      hyp_lt_mc_motherid  ->push_back(mus_mc_motherid  ->at(tight_index)  );
      hyp_lt_index        ->push_back(tight_index                         );
      hyp_lt_id           ->push_back(-13*(mus_charge   ->at(tight_index)) );
      hyp_lt_d0           ->push_back(mus_d0           ->at(tight_index)  );
      hyp_lt_z0           ->push_back(mus_z0           ->at(tight_index)  );
      hyp_lt_d0corr       ->push_back(mus_d0corr       ->at(tight_index)  );
      hyp_lt_z0corr       ->push_back(mus_z0corr       ->at(tight_index)  );
      hyp_lt_vertexphi    ->push_back(mus_vertexphi    ->at(tight_index)  );
      hyp_lt_chi2         ->push_back(mus_chi2         ->at(tight_index)  );
      hyp_lt_ndof         ->push_back(mus_ndof         ->at(tight_index)  );
      hyp_lt_d0Err        ->push_back(mus_d0Err        ->at(tight_index)  );
      hyp_lt_z0Err        ->push_back(mus_z0Err        ->at(tight_index)  );
      hyp_lt_ptErr        ->push_back(mus_ptErr        ->at(tight_index)  );
      hyp_lt_etaErr       ->push_back(mus_etaErr       ->at(tight_index)  );
      hyp_lt_phiErr       ->push_back(mus_phiErr       ->at(tight_index)  );
      hyp_lt_outerPhi     ->push_back(mus_outerPhi     ->at(tight_index)  );
      hyp_lt_outerEta     ->push_back(mus_outerEta     ->at(tight_index)  );
      hyp_lt_iso          ->push_back(mus_iso          ->at(tight_index)  );
      hyp_lt_tkIso        ->push_back(mus_iso03_sumPt  ->at(tight_index)  );
      hyp_lt_p4           ->push_back(mus_p4           ->at(tight_index)  );
      hyp_lt_trk_p4       ->push_back(mus_trk_p4       ->at(tight_index)  );
      hyp_lt_mc_p4        ->push_back(mus_mc_p4        ->at(tight_index)  );
      
      hyp_ll_validHits    ->push_back(mus_validHits    ->at(loose_index)  );
      hyp_ll_lostHits     ->push_back(mus_lostHits     ->at(loose_index)  );
      hyp_ll_mc_id        ->push_back(mus_mc_id        ->at(loose_index)  );
      hyp_ll_charge       ->push_back(mus_charge       ->at(loose_index)  );
      hyp_ll_mc_motherid  ->push_back(mus_mc_motherid  ->at(loose_index)  );
      hyp_ll_index        ->push_back(loose_index                         );
      hyp_ll_id           ->push_back(-13*(mus_charge   ->at(loose_index)) );
      hyp_ll_d0           ->push_back(mus_d0           ->at(loose_index)  );
      hyp_ll_z0           ->push_back(mus_z0           ->at(loose_index)  );
      hyp_ll_d0corr       ->push_back(mus_d0corr       ->at(loose_index)  );
      hyp_ll_z0corr       ->push_back(mus_z0corr       ->at(loose_index)  );
      hyp_ll_vertexphi    ->push_back(mus_vertexphi    ->at(loose_index)  );
      hyp_ll_chi2         ->push_back(mus_chi2         ->at(loose_index)  );
      hyp_ll_ndof         ->push_back(mus_ndof         ->at(loose_index)  );
      hyp_ll_d0Err        ->push_back(mus_d0Err        ->at(loose_index)  );
      hyp_ll_z0Err        ->push_back(mus_z0Err        ->at(loose_index)  );
      hyp_ll_ptErr        ->push_back(mus_ptErr        ->at(loose_index)  );
      hyp_ll_etaErr       ->push_back(mus_etaErr       ->at(loose_index)  );
      hyp_ll_phiErr       ->push_back(mus_phiErr       ->at(loose_index)  );
      hyp_ll_outerPhi     ->push_back(mus_outerPhi     ->at(loose_index)  );
      hyp_ll_outerEta     ->push_back(mus_outerEta     ->at(loose_index)  );
      hyp_ll_iso          ->push_back(mus_iso          ->at(loose_index)  );
      hyp_ll_tkIso        ->push_back(mus_iso03_sumPt  ->at(loose_index)  );
      hyp_ll_p4           ->push_back(mus_p4           ->at(loose_index)  );
      hyp_ll_trk_p4       ->push_back(mus_trk_p4       ->at(loose_index)  );
      hyp_ll_mc_p4        ->push_back(mus_mc_p4        ->at(loose_index)  );

      hyp_met             ->push_back(hypmet                                 );
      hyp_metPhi          ->push_back(hypmetPhi                              );
      hyp_metCaloExp      ->push_back(hypmet_MIP                             );
      hyp_metPhiCaloExp   ->push_back(hypmetPhi_MIP                          );
      hyp_metCone         ->push_back(hypmet_S9                              );
      hyp_metPhiCone      ->push_back(hypmetPhi_S9                           );
      hyp_metNoCalo       ->push_back(hypmet_nocalo                          );
      hyp_metPhiNoCalo    ->push_back(hypmetPhi_nocalo                       );
      hyp_metAll          ->push_back(metAll                                 );
      hyp_metPhiAll       ->push_back(metPhiAll                              );
      hyp_metAllCaloExp   ->push_back(metAllCaloExp                          );
      hyp_metPhiAllCaloExp->push_back(metPhiAllCaloExp                       );
      hyp_metJes5         ->push_back(metJes5                                );
      hyp_metPhiJes5      ->push_back(metPhiJes5                             );
      hyp_metJes10        ->push_back(metJes10                               );
      hyp_metPhiJes10     ->push_back(metPhiJes10                            );
      hyp_metJes15        ->push_back(metJes15                               );
      hyp_metPhiJes15     ->push_back(metPhiJes15                            );
      hyp_metJes30        ->push_back(metJes30                               );
      hyp_metPhiJes30     ->push_back(metPhiJes30                            );
      hyp_metJes50        ->push_back(metJes50                               );
      hyp_metPhiJes50     ->push_back(metPhiJes50                            );
      //   hyp_metEMF5         ->push_back(metEMF5                                );
      //       hyp_metPhiEMF5      ->push_back(metPhiEMF5                             );
      //       hyp_metEMF10        ->push_back(metEMF10                               );
      //       hyp_metPhiEMF10     ->push_back(metPhiEMF10                            );
      //       hyp_metEMF15        ->push_back(metEMF15                               );
      //       hyp_metPhiEMF15     ->push_back(metPhiEMF15                            );
      //       hyp_metEMF30        ->push_back(metEMF30                               );
      //       hyp_metPhiEMF30     ->push_back(metPhiEMF30                            );
      //       hyp_metEMF50        ->push_back(metEMF50                               );
      //       hyp_metPhiEMF50     ->push_back(metPhiEMF50                            );
      hyp_metDPhiJet10    ->push_back(METUtilities::metObjDPhi(*jets_p4,
							       hypmetPhi,
							       10.)         );
      hyp_metDPhiJet15    ->push_back(METUtilities::metObjDPhi(*jets_p4,
							       hypmetPhi,
							       15.)         ); 
      hyp_metDPhiJet20    ->push_back(METUtilities::metObjDPhi(*jets_p4,
							       hypmetPhi,
							       20.)         );
      hyp_metDPhiTrk10    ->push_back(METUtilities::metObjDPhi(*trks_p4,
							       hypmetPhi,
							       10.)          );
      hyp_metDPhiTrk25    ->push_back(METUtilities::metObjDPhi(*trks_p4,
							       hypmetPhi,
							       25.)         );
      hyp_metDPhiTrk50    ->push_back(METUtilities::metObjDPhi(*trks_p4,
							       hypmetPhi,
							       50.)         );
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
      vector<int>            temp_jets_mc_id;
      vector<float>          temp_jets_emFrac;
      vector<float>          temp_jets_chFrac;
      vector<float>          temp_jets_mc_emEnergy;
      vector<float>          temp_jets_mc_hadEnergy;
      vector<float>          temp_jets_mc_invEnergy;
      vector<float>          temp_jets_mc_otherEnergy;
      vector<float>          temp_jets_cor;
      vector<float>          temp_jets_EMFcor;
                                                   
      vector<LorentzVector>  temp_jets_p4;
      vector<LorentzVector>  temp_jets_mc_p4;
      vector<LorentzVector>  temp_jets_mc_gp_p4;

      vector<int>            temp_other_jets_mc_id;
      vector<float>          temp_other_jets_emFrac;
      vector<float>          temp_other_jets_chFrac;
      vector<float>          temp_other_jets_mc_emEnergy;
      vector<float>          temp_other_jets_mc_hadEnergy;
      vector<float>          temp_other_jets_mc_invEnergy;
      vector<float>          temp_other_jets_mc_otherEnergy;
      vector<float>          temp_other_jets_cor;
      vector<float>          temp_other_jets_EMFcor;
      
      vector<LorentzVector>  temp_other_jets_p4;
      vector<LorentzVector>  temp_other_jets_mc_p4;
      vector<LorentzVector>  temp_other_jets_mc_gp_p4;

      vector<int>            temp_jets_pat_genParton_id;
      vector<int>            temp_jets_pat_genPartonMother_id;
      vector<int>            temp_jets_pat_partonFlavour;
      vector<float>          temp_jets_pat_noCorrF;                  
      vector<float>          temp_jets_pat_udsCorrF;                 
      vector<float>          temp_jets_pat_gluCorrF;                 
      vector<float>          temp_jets_pat_cCorrF;                   
      vector<float>          temp_jets_pat_bCorrF;                   
      vector<float>          temp_jets_pat_jetCharge;   
      vector<LorentzVector>  temp_jets_pat_genParton_p4;
      vector<LorentzVector>  temp_jets_pat_genPartonMother_p4;
      vector<LorentzVector>  temp_jets_pat_genJet_p4;
      vector<LorentzVector>  temp_jets_pat_jet_p4;
      
      vector<int>            temp_other_jets_pat_genParton_id;
      vector<int>            temp_other_jets_pat_genPartonMother_id;
      vector<int>            temp_other_jets_pat_partonFlavour;
      vector<float>          temp_other_jets_pat_noCorrF;                  
      vector<float>          temp_other_jets_pat_udsCorrF;                 
      vector<float>          temp_other_jets_pat_gluCorrF;                 
      vector<float>          temp_other_jets_pat_cCorrF;                   
      vector<float>          temp_other_jets_pat_bCorrF;                   
      vector<float>          temp_other_jets_pat_jetCharge;     
      vector<LorentzVector>  temp_other_jets_pat_genParton_p4;
      vector<LorentzVector>  temp_other_jets_pat_genPartonMother_p4;
      vector<LorentzVector>  temp_other_jets_pat_genJet_p4;
      vector<LorentzVector>  temp_other_jets_pat_jet_p4;

      //these are for the MET correction later
      vector<LorentzVector> jets_noel_p4;
      vector<LorentzVector> jets_noelEMF_p4;
      vector<float> jets_noel_jescor;
      //vector<float> jets_noel_EMFcor;
      
      for(unsigned int i = 0; i<jets_p4->size(); i++) {

	// we don't want jets that overlap with electrons
	if(!testJetForElectrons(jets_p4->at(i), els_p4->at(loose_index))) continue;
	if(!testJetForElectrons(jets_p4->at(i), els_p4->at(tight_index))) continue;
	
	float px = (jets_p4->at(i).px())*jets_pat_noCorrF_h->at(i);
	float py = (jets_p4->at(i).py())*jets_pat_noCorrF_h->at(i);
	float pz = (jets_p4->at(i).pz())*jets_pat_noCorrF_h->at(i);
	float E  = (jets_p4->at(i).E())*jets_pat_noCorrF_h->at(i);
	LorentzVector temp(px, py, pz, E);
	if(usingPATJets) {
	  jets_noel_p4        .push_back(temp); //uncorrected p4
	  jets_noel_jescor    .push_back(jets_cor->at(i)/(jets_pat_noCorrF_h->at(i)));
	} else {
	  jets_noel_p4        .push_back(jets_p4    ->at(i));
	  jets_noel_jescor    .push_back(jets_cor   ->at(i));
	}
	
	px = (jets_EMFcor->at(i))*(jets_p4->at(i).px());
	py = (jets_EMFcor->at(i))*(jets_p4->at(i).py());
	pz = (jets_EMFcor->at(i))*(jets_p4->at(i).pz());
	E =  (jets_EMFcor->at(i))*(jets_p4->at(i).E());
	temp = LorentzVector(px, py, pz, E);
	//jets_noelEMF_p4.push_back(temp);
	//jets_noel_EMFcor.push_back(jets_EMFcor->at(i));
	  

	double jet_eta = jets_p4->at(i).eta();
	double jet_pt = jets_p4->at(i).Pt();
	double jet_ptcut;
	
	if(usingPATJets) {
	  jet_ptcut = hypJetMinPtCut/(jets_pat_noCorrF_h->at(i));
	} else jet_ptcut = hypJetMinPtCut;
	
	if( hypJetMinEtaCut < jet_eta && 
	    jet_eta < hypJetMaxEtaCut && 
	    jet_pt  > jet_ptcut ) { //hyp jets
	  
	  temp_jets_mc_id                  .push_back(jets_mc_id           ->at(i));
	  temp_jets_emFrac                 .push_back(jets_emFrac          ->at(i));
	  temp_jets_chFrac                 .push_back(jets_chFrac          ->at(i));
	  
	  temp_jets_mc_emEnergy            .push_back(jets_mc_emEnergy     ->at(i));
	  temp_jets_mc_hadEnergy           .push_back(jets_mc_hadEnergy    ->at(i));
	  temp_jets_mc_invEnergy           .push_back(jets_mc_invEnergy    ->at(i));
	  temp_jets_mc_otherEnergy         .push_back(jets_mc_otherEnergy  ->at(i));  
	  temp_jets_cor                    .push_back(jets_cor             ->at(i));
	  temp_jets_EMFcor                 .push_back(jets_EMFcor          ->at(i));
	  
	  temp_jets_p4                     .push_back(jets_p4              ->at(i));
	  temp_jets_mc_p4                  .push_back(jets_mc_p4           ->at(i));
	  temp_jets_mc_gp_p4               .push_back(jets_mc_gp_p4        ->at(i));

	   
	  if(usingPATJets) {
	    temp_jets_pat_genParton_id      .push_back(jets_pat_genParton_id_h       ->at(i));
	    temp_jets_pat_genPartonMother_id.push_back(jets_pat_genPartonMother_id_h ->at(i));
	    temp_jets_pat_partonFlavour     .push_back(jets_pat_partonFlavour_h      ->at(i));
	    temp_jets_pat_genParton_p4      .push_back(jets_pat_genParton_p4_h       ->at(i));
	    temp_jets_pat_genPartonMother_p4.push_back(jets_pat_genPartonMother_p4_h ->at(i));
	    temp_jets_pat_genJet_p4         .push_back(jets_pat_genJet_p4_h          ->at(i));
	    temp_jets_pat_jet_p4            .push_back(jets_pat_jet_p4_h             ->at(i));
	    temp_jets_pat_noCorrF           .push_back(jets_pat_noCorrF_h            ->at(i));                  
	    temp_jets_pat_udsCorrF          .push_back(jets_pat_udsCorrF_h           ->at(i));                 
	    temp_jets_pat_gluCorrF          .push_back(jets_pat_gluCorrF_h           ->at(i));                 
	    temp_jets_pat_cCorrF            .push_back(jets_pat_cCorrF_h             ->at(i));                   
	    temp_jets_pat_bCorrF            .push_back(jets_pat_bCorrF_h             ->at(i));                   
	    temp_jets_pat_jetCharge        .push_back(jets_pat_jetCharge_h            ->at(i));               
	  }
	  
      } else {
	temp_other_jets_mc_id            .push_back(jets_mc_id           ->at(i));
	temp_other_jets_emFrac           .push_back(jets_emFrac          ->at(i));
	temp_other_jets_chFrac           .push_back(jets_chFrac          ->at(i));
	  
	temp_other_jets_mc_emEnergy      .push_back(jets_mc_emEnergy     ->at(i));
	temp_other_jets_mc_hadEnergy     .push_back(jets_mc_hadEnergy    ->at(i));
	temp_other_jets_mc_invEnergy     .push_back(jets_mc_invEnergy    ->at(i));
	temp_other_jets_mc_otherEnergy   .push_back(jets_mc_otherEnergy  ->at(i));  
	temp_other_jets_cor              .push_back(jets_cor             ->at(i));
	temp_other_jets_EMFcor           .push_back(jets_EMFcor          ->at(i));
	  
	temp_other_jets_p4               .push_back(jets_p4              ->at(i));
	temp_other_jets_mc_p4            .push_back(jets_mc_p4           ->at(i));
	temp_other_jets_mc_gp_p4         .push_back(jets_mc_gp_p4        ->at(i));
  
	if(usingPATJets) {
	  temp_other_jets_pat_genParton_id      .push_back(jets_pat_genParton_id_h       ->at(i));
	  temp_other_jets_pat_genPartonMother_id.push_back(jets_pat_genPartonMother_id_h ->at(i));
	  temp_other_jets_pat_partonFlavour     .push_back(jets_pat_partonFlavour_h      ->at(i));
	  temp_other_jets_pat_genParton_p4      .push_back(jets_pat_genParton_p4_h       ->at(i));
	  temp_other_jets_pat_genPartonMother_p4.push_back(jets_pat_genPartonMother_p4_h ->at(i));
	  temp_other_jets_pat_genJet_p4         .push_back(jets_pat_genJet_p4_h          ->at(i));
	  temp_other_jets_pat_jet_p4            .push_back(jets_pat_jet_p4_h             ->at(i));
	  temp_other_jets_pat_noCorrF           .push_back(jets_pat_noCorrF_h            ->at(i));                  
	  temp_other_jets_pat_udsCorrF          .push_back(jets_pat_udsCorrF_h           ->at(i));                 
	  temp_other_jets_pat_gluCorrF          .push_back(jets_pat_gluCorrF_h           ->at(i));                 
	  temp_other_jets_pat_cCorrF            .push_back(jets_pat_cCorrF_h             ->at(i));                   
	  temp_other_jets_pat_bCorrF            .push_back(jets_pat_bCorrF_h             ->at(i));                   
	  temp_other_jets_pat_jetCharge         .push_back(jets_pat_jetCharge_h            ->at(i));               
	}
	  
      }
    }

    //push these into the hyp_jets and hyp_other_jets vars
    hyp_jets_mc_id                 ->push_back(temp_jets_mc_id              );
    hyp_jets_emFrac                ->push_back(temp_jets_emFrac             );
    hyp_jets_chFrac                ->push_back(temp_jets_chFrac             );
    hyp_jets_mc_emEnergy           ->push_back(temp_jets_mc_emEnergy        );
    hyp_jets_mc_hadEnergy          ->push_back(temp_jets_mc_hadEnergy       );
    hyp_jets_mc_invEnergy          ->push_back(temp_jets_mc_invEnergy       );
    hyp_jets_mc_otherEnergy        ->push_back(temp_jets_mc_otherEnergy     );
    hyp_jets_cor                   ->push_back(temp_jets_cor                );
    hyp_jets_EMFcor                ->push_back(temp_jets_EMFcor             );
      
    hyp_jets_p4                    ->push_back(temp_jets_p4                  );
    hyp_jets_mc_p4                 ->push_back(temp_jets_mc_p4               );
    hyp_jets_mc_gp_p4              ->push_back(temp_jets_mc_gp_p4            );
      
    hyp_other_jets_mc_id           ->push_back(temp_other_jets_mc_id         );
    hyp_other_jets_emFrac          ->push_back(temp_other_jets_emFrac        );
    hyp_other_jets_chFrac          ->push_back(temp_other_jets_chFrac        );
    hyp_other_jets_mc_emEnergy     ->push_back(temp_other_jets_mc_emEnergy   );
    hyp_other_jets_mc_hadEnergy    ->push_back(temp_other_jets_mc_hadEnergy  );
    hyp_other_jets_mc_invEnergy    ->push_back(temp_other_jets_mc_invEnergy  );
    hyp_other_jets_mc_otherEnergy  ->push_back(temp_other_jets_mc_otherEnergy);
    hyp_other_jets_cor             ->push_back(temp_other_jets_cor           );
    hyp_other_jets_EMFcor          ->push_back(temp_other_jets_EMFcor        );
      
    hyp_other_jets_p4              ->push_back(temp_other_jets_p4            );
    hyp_other_jets_mc_p4           ->push_back(temp_other_jets_mc_p4         );
    hyp_other_jets_mc_gp_p4        ->push_back(temp_other_jets_mc_gp_p4      );

    if(usingPATJets) {
      hyp_jets_pat_genParton_id            ->push_back(temp_jets_pat_genParton_id             );
      hyp_jets_pat_genPartonMother_id      ->push_back(temp_jets_pat_genPartonMother_id       );
      hyp_jets_pat_partonFlavour           ->push_back(temp_jets_pat_partonFlavour            );
      hyp_jets_pat_genParton_p4            ->push_back(temp_jets_pat_genParton_p4             );
      hyp_jets_pat_genPartonMother_p4      ->push_back(temp_jets_pat_genPartonMother_p4       );
      hyp_jets_pat_genJet_p4               ->push_back(temp_jets_pat_genJet_p4                );
      hyp_jets_pat_noCorrF                 ->push_back(temp_jets_pat_noCorrF                );                  
      hyp_jets_pat_udsCorrF                ->push_back(temp_jets_pat_udsCorrF               );                 
      hyp_jets_pat_gluCorrF                ->push_back(temp_jets_pat_gluCorrF               );                 
      hyp_jets_pat_cCorrF                  ->push_back(temp_jets_pat_cCorrF                 );                   
      hyp_jets_pat_bCorrF                  ->push_back(temp_jets_pat_bCorrF                 );                   
      hyp_jets_pat_jetCharge              ->push_back(temp_jets_pat_jetCharge                );               

      hyp_jets_pat_jet_p4                  ->push_back(temp_jets_pat_jet_p4                   );
      hyp_other_jets_pat_genParton_id      ->push_back(temp_other_jets_pat_genParton_id       );
      hyp_other_jets_pat_genPartonMother_id->push_back(temp_other_jets_pat_genPartonMother_id );
      hyp_other_jets_pat_partonFlavour     ->push_back(temp_other_jets_pat_partonFlavour      );
      hyp_other_jets_pat_genParton_p4      ->push_back(temp_other_jets_pat_genParton_p4       );
      hyp_other_jets_pat_genPartonMother_p4->push_back(temp_other_jets_pat_genPartonMother_p4 );
      hyp_other_jets_pat_genJet_p4         ->push_back(temp_other_jets_pat_genJet_p4          );
      hyp_other_jets_pat_jet_p4            ->push_back(temp_other_jets_pat_jet_p4             );
      hyp_other_jets_pat_noCorrF           ->push_back(temp_other_jets_pat_noCorrF          );                  
      hyp_other_jets_pat_udsCorrF          ->push_back(temp_other_jets_pat_udsCorrF         );                 
      hyp_other_jets_pat_gluCorrF          ->push_back(temp_other_jets_pat_gluCorrF         );                 
      hyp_other_jets_pat_cCorrF            ->push_back(temp_other_jets_pat_cCorrF           );                   
      hyp_other_jets_pat_bCorrF            ->push_back(temp_other_jets_pat_bCorrF           );                   
      hyp_other_jets_pat_jetCharge        ->push_back(temp_other_jets_pat_jetCharge          );               
    }

	
    double metJes5 = *evt_met;
    double metPhiJes5 = *evt_metphi;
    double metJes10 = *evt_met;
    double metPhiJes10 = *evt_metphi;
    double metJes15 = *evt_met;
    double metPhiJes15 = *evt_metphi;
    double metJes30 = *evt_met;
    double metPhiJes30 = *evt_metphi;
    double metJes50 = *evt_met;
    double metPhiJes50 = *evt_metphi;
       
    //       double metEMF5 = *evt_met;
    //       double metPhiEMF5 = *evt_metphi;
    //       double metEMF10 = *evt_met;
    //       double metPhiEMF10 = *evt_metphi;
    //       double metEMF15 = *evt_met;
    //       double metPhiEMF15 = *evt_metphi;
    //       double metEMF30 = *evt_met;
    //       double metPhiEMF30 = *evt_metphi;
    //       double metEMF50 = *evt_met;
    //       double metPhiEMF50 = *evt_metphi;
       
     
  
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes5, metPhiJes5, 5);
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes10, metPhiJes10, 10);
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes15, metPhiJes15, 15);
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes30, metPhiJes30, 30);
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes50, metPhiJes50, 50);
      
             
    //       METUtilities::correctedJetMET(jets_noelEMF_p4, jets_noel_EMFcor,
    // 				    metEMF5, metPhiEMF5, 5);
    //       METUtilities::correctedJetMET(jets_noelEMF_p4, jets_noel_EMFcor,
    // 				    metEMF10, metPhiEMF10, 10);
    //       METUtilities::correctedJetMET(jets_noelEMF_p4, jets_noel_EMFcor,
    // 				    metEMF15, metPhiEMF15, 15);
    //       METUtilities::correctedJetMET(jets_noelEMF_p4, jets_noel_EMFcor,
    // 				    metEMF30, metPhiEMF30, 30);
    //       METUtilities::correctedJetMET(jets_noelEMF_p4, jets_noel_EMFcor,
    // 				    metEMF50, metPhiEMF50, 50);



    hyp_type          ->push_back(3);
    hyp_njets         ->push_back(temp_jets_p4.size()                   );     
    hyp_nojets        ->push_back(temp_other_jets_p4.size()             );     
    hyp_p4            ->push_back(els_p4->at(tight_index)+
				  els_p4->at(loose_index)               );
      
    hyp_lt_validHits    ->push_back(els_validHits    ->at(tight_index)  );
    hyp_lt_lostHits     ->push_back(els_lostHits     ->at(tight_index)  );
    hyp_lt_mc_id        ->push_back(els_mc_id        ->at(tight_index)  );
    hyp_lt_charge       ->push_back(els_charge       ->at(tight_index)  );
    hyp_lt_mc_motherid  ->push_back(els_mc_motherid  ->at(tight_index)  );
    hyp_lt_index        ->push_back(tight_index                         );
    hyp_lt_id           ->push_back(-11*(els_charge   ->at(tight_index)));
    hyp_lt_d0           ->push_back(els_d0           ->at(tight_index)  );
    hyp_lt_z0           ->push_back(els_z0           ->at(tight_index)  );
    hyp_lt_d0corr       ->push_back(els_d0corr       ->at(tight_index)  );
    hyp_lt_z0corr       ->push_back(els_z0corr       ->at(tight_index)  );
    hyp_lt_vertexphi    ->push_back(els_vertexphi    ->at(tight_index)  );
    hyp_lt_chi2         ->push_back(els_chi2         ->at(tight_index)  );
    hyp_lt_ndof         ->push_back(els_ndof         ->at(tight_index)  );
    hyp_lt_d0Err        ->push_back(els_d0Err        ->at(tight_index)  );
    hyp_lt_z0Err        ->push_back(els_z0Err        ->at(tight_index)  );
    hyp_lt_ptErr        ->push_back(els_ptErr        ->at(tight_index)  );
    hyp_lt_etaErr       ->push_back(els_etaErr       ->at(tight_index)  );
    hyp_lt_phiErr       ->push_back(els_phiErr       ->at(tight_index)  );
    hyp_lt_outerPhi     ->push_back(els_outerPhi     ->at(tight_index)  );
    hyp_lt_outerEta     ->push_back(els_outerEta     ->at(tight_index)  );
    hyp_lt_iso          ->push_back(els_tkIso        ->at(tight_index)  );
    hyp_lt_tkIso        ->push_back(els_tkIso        ->at(tight_index)  );
    hyp_lt_p4           ->push_back(els_p4           ->at(tight_index)  );
    hyp_lt_trk_p4       ->push_back(els_trk_p4       ->at(tight_index)  );
    hyp_lt_mc_p4        ->push_back(els_mc_p4        ->at(tight_index)  );
            
    hyp_ll_validHits    ->push_back(els_validHits    ->at(loose_index)  );
    hyp_ll_lostHits     ->push_back(els_lostHits     ->at(loose_index)  );
    hyp_ll_mc_id        ->push_back(els_mc_id        ->at(loose_index)  );
    hyp_ll_charge       ->push_back(els_charge       ->at(loose_index)  );
    hyp_ll_mc_motherid  ->push_back(els_mc_motherid  ->at(loose_index)  );
    hyp_ll_index        ->push_back(loose_index                         );
    hyp_ll_id           ->push_back(-11*(els_charge   ->at(loose_index)) );
    hyp_ll_d0           ->push_back(els_d0           ->at(loose_index)  );
    hyp_ll_z0           ->push_back(els_z0           ->at(loose_index)  );
    hyp_ll_d0corr       ->push_back(els_d0corr       ->at(loose_index)  );
    hyp_ll_z0corr       ->push_back(els_z0corr       ->at(loose_index)  );
    hyp_ll_vertexphi    ->push_back(els_vertexphi    ->at(loose_index)  );
    hyp_ll_chi2         ->push_back(els_chi2         ->at(loose_index)  );
    hyp_ll_ndof         ->push_back(els_ndof         ->at(loose_index)  );
    hyp_ll_d0Err        ->push_back(els_d0Err        ->at(loose_index)  );
    hyp_ll_z0Err        ->push_back(els_z0Err        ->at(loose_index)  );
    hyp_ll_ptErr        ->push_back(els_ptErr        ->at(loose_index)  );
    hyp_ll_etaErr       ->push_back(els_etaErr       ->at(loose_index)  );
    hyp_ll_phiErr       ->push_back(els_phiErr       ->at(loose_index)  );
    hyp_ll_outerPhi     ->push_back(els_outerPhi     ->at(loose_index)  );
    hyp_ll_outerEta     ->push_back(els_outerEta     ->at(loose_index)  );
    hyp_ll_iso          ->push_back(els_tkIso        ->at(loose_index)  );
    hyp_ll_tkIso        ->push_back(els_tkIso        ->at(loose_index)  );
    hyp_ll_p4           ->push_back(els_p4           ->at(loose_index)  );
    hyp_ll_trk_p4       ->push_back(els_trk_p4       ->at(loose_index)  );
    hyp_ll_mc_p4        ->push_back(els_mc_p4        ->at(loose_index)  );

    hyp_met             ->push_back(*evt_met                            );
    hyp_metPhi          ->push_back(*evt_metphi                         );
    hyp_metCaloExp      ->push_back(*evt_met                            );
    hyp_metPhiCaloExp   ->push_back(*evt_metphi                         );
    hyp_metCone         ->push_back(*evt_met                            );
    hyp_metPhiCone      ->push_back(*evt_metphi                         );
    hyp_metNoCalo       ->push_back(*evt_met                            );
    hyp_metPhiNoCalo    ->push_back(*evt_metphi                         );
    hyp_metAll          ->push_back(metAll                              );
    hyp_metPhiAll       ->push_back(metPhiAll                           );
    hyp_metAllCaloExp   ->push_back(metAllCaloExp                       );
    hyp_metPhiAllCaloExp->push_back(metPhiAllCaloExp                    );
    hyp_metJes5         ->push_back(metJes5                             );
    hyp_metPhiJes5      ->push_back(metPhiJes5                          );
    hyp_metJes10        ->push_back(metJes10                            );
    hyp_metPhiJes10     ->push_back(metPhiJes10                         );
    hyp_metJes15        ->push_back(metJes15                            );
    hyp_metPhiJes15     ->push_back(metPhiJes15                         );
    hyp_metJes30        ->push_back(metJes30                            );
    hyp_metPhiJes30     ->push_back(metPhiJes30                         );
    hyp_metJes50        ->push_back(metJes50                            );
    hyp_metPhiJes50     ->push_back(metPhiJes50                         );
    //   hyp_metEMF5         ->push_back(metEMF5                             );
    //       hyp_metPhiEMF5      ->push_back(metPhiEMF5                          );
    //       hyp_metEMF10        ->push_back(metEMF10                            );
    //       hyp_metPhiEMF10     ->push_back(metPhiEMF10                         );
    //       hyp_metEMF15        ->push_back(metEMF15                            );
    //       hyp_metPhiEMF15     ->push_back(metPhiEMF15                         );
    //       hyp_metEMF30        ->push_back(metEMF30                            );
    //       hyp_metPhiEMF30     ->push_back(metPhiEMF30                         );
    //       hyp_metEMF50        ->push_back(metEMF50                            );
    //       hyp_metPhiEMF50     ->push_back(metPhiEMF50                         );
    hyp_metDPhiJet10    ->push_back(METUtilities::metObjDPhi(*jets_p4,
							     *evt_metphi,
							     10.)       );
    hyp_metDPhiJet15    ->push_back(METUtilities::metObjDPhi(*jets_p4,
							     *evt_metphi,
							     15.)       ); 
    hyp_metDPhiJet20    ->push_back(METUtilities::metObjDPhi(*jets_p4,
							     *evt_metphi,
							     20.)       );
    hyp_metDPhiTrk10    ->push_back(METUtilities::metObjDPhi(*trks_p4,
							     *evt_metphi,
							     10.)       );
    hyp_metDPhiTrk25    ->push_back(METUtilities::metObjDPhi(*trks_p4,
							     *evt_metphi,
							     25.)       );
    hyp_metDPhiTrk50    ->push_back(METUtilities::metObjDPhi(*trks_p4,
							     *evt_metphi,
							     50.)       );

  }
}  



  
/*------------------------------------------------------------
  The EMu, MuE cases
  To avoid double counting, only make MuE if Mu is tight and E is loose
*/

for(unsigned int els_index = 0; els_index < nels; els_index++) {
  for(unsigned int mus_index = 0; mus_index < nmus; mus_index++) {

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
      
    //fill the Jet vars
    vector<int>            temp_jets_mc_id;
    vector<float>          temp_jets_emFrac;
    vector<float>          temp_jets_chFrac;
    vector<float>          temp_jets_mc_emEnergy;
    vector<float>          temp_jets_mc_hadEnergy;
    vector<float>          temp_jets_mc_invEnergy;
    vector<float>          temp_jets_mc_otherEnergy;
    vector<float>          temp_jets_cor;
    vector<float>          temp_jets_EMFcor;
                                                   
    vector<LorentzVector>  temp_jets_p4;
    vector<LorentzVector>  temp_jets_mc_p4;
    vector<LorentzVector>  temp_jets_mc_gp_p4;

    vector<int>            temp_other_jets_mc_id;
    vector<float>          temp_other_jets_emFrac;
    vector<float>          temp_other_jets_chFrac;
    vector<float>          temp_other_jets_mc_emEnergy;
    vector<float>          temp_other_jets_mc_hadEnergy;
    vector<float>          temp_other_jets_mc_invEnergy;
    vector<float>          temp_other_jets_mc_otherEnergy;
    vector<float>          temp_other_jets_cor;
    vector<float>          temp_other_jets_EMFcor;
      
    vector<LorentzVector>  temp_other_jets_p4;
    vector<LorentzVector>  temp_other_jets_mc_p4;
    vector<LorentzVector>  temp_other_jets_mc_gp_p4;
      
    vector<int>            temp_jets_pat_genParton_id;
    vector<int>            temp_jets_pat_genPartonMother_id;
    vector<int>            temp_jets_pat_partonFlavour;
    vector<float>          temp_jets_pat_noCorrF;                  
    vector<float>          temp_jets_pat_udsCorrF;                 
    vector<float>          temp_jets_pat_gluCorrF;                 
    vector<float>          temp_jets_pat_cCorrF;                   
    vector<float>          temp_jets_pat_bCorrF;                   
    vector<float>          temp_jets_pat_jetCharge;   
    vector<LorentzVector>  temp_jets_pat_genParton_p4;
    vector<LorentzVector>  temp_jets_pat_genPartonMother_p4;
    vector<LorentzVector>  temp_jets_pat_genJet_p4;
    vector<LorentzVector>  temp_jets_pat_jet_p4;
      
    vector<int>            temp_other_jets_pat_genParton_id;
    vector<int>            temp_other_jets_pat_genPartonMother_id;
    vector<int>            temp_other_jets_pat_partonFlavour;
    vector<float>          temp_other_jets_pat_noCorrF;                  
    vector<float>          temp_other_jets_pat_udsCorrF;                 
    vector<float>          temp_other_jets_pat_gluCorrF;                 
    vector<float>          temp_other_jets_pat_cCorrF;                   
    vector<float>          temp_other_jets_pat_bCorrF;                   
    vector<float>          temp_other_jets_pat_jetCharge;    
    vector<LorentzVector>  temp_other_jets_pat_genParton_p4;
    vector<LorentzVector>  temp_other_jets_pat_genPartonMother_p4;
    vector<LorentzVector>  temp_other_jets_pat_genJet_p4;
    vector<LorentzVector>  temp_other_jets_pat_jet_p4;
      
    //these are for the MET correction later
    vector<LorentzVector> jets_noel_p4;
    vector<float> jets_noel_jescor;
    //vector<float> jets_noel_EMFcor;
    //for the jet corrected MET, we need all jets that do not
    //overlap with the hyp electrons
       
    for(unsigned int i = 0; i<jets_p4->size(); i++) {
		
      //we don't any jets that overlap with an electron
      if(!testJetForElectrons(jets_p4->at(i), els_p4->at(els_index))) continue;
	
      if(usingPATJets) {
	float px = (jets_p4->at(i).px())*jets_pat_noCorrF_h->at(i);
	float py = (jets_p4->at(i).py())*jets_pat_noCorrF_h->at(i);
	float pz = (jets_p4->at(i).pz())*jets_pat_noCorrF_h->at(i);
	float E  = (jets_p4->at(i).E())*jets_pat_noCorrF_h->at(i);
	LorentzVector temp(px, py, pz, E);
	jets_noel_p4        .push_back(temp);
	jets_noel_jescor    .push_back(jets_cor->at(i)/(jets_pat_noCorrF_h->at(i)));
      } else {
	jets_noel_p4        .push_back(jets_p4    ->at(i));
	jets_noel_jescor    .push_back(jets_cor   ->at(i));
      }

      //jets_noel_EMFcor    .push_back(jets_EMFcor->at(i));  

      double jet_eta = jets_p4->at(i).eta();
      double jet_pt  = jets_p4->at(i).Pt();
      double jet_ptcut;

      if(usingPATJets) {
	jet_ptcut = hypJetMinPtCut/(jets_pat_noCorrF_h->at(i));
      } else jet_ptcut =  hypJetMinPtCut;
	
      if( hypJetMinEtaCut < jet_eta && 
	  jet_eta < hypJetMaxEtaCut && 
	  jet_pt  > jet_ptcut) { 
	  
	temp_jets_mc_id                  .push_back(jets_mc_id           ->at(i));
	temp_jets_emFrac                 .push_back(jets_emFrac          ->at(i));
	temp_jets_chFrac                 .push_back(jets_chFrac          ->at(i));
	  
	temp_jets_mc_emEnergy            .push_back(jets_mc_emEnergy     ->at(i));
	temp_jets_mc_hadEnergy           .push_back(jets_mc_hadEnergy    ->at(i));
	temp_jets_mc_invEnergy           .push_back(jets_mc_invEnergy    ->at(i));
	temp_jets_mc_otherEnergy         .push_back(jets_mc_otherEnergy  ->at(i));  
	temp_jets_cor                    .push_back(jets_cor             ->at(i));
	temp_jets_EMFcor                 .push_back(jets_EMFcor          ->at(i));
	  
	temp_jets_p4                     .push_back(jets_p4              ->at(i));
	temp_jets_mc_p4                  .push_back(jets_mc_p4           ->at(i));
	temp_jets_mc_gp_p4               .push_back(jets_mc_gp_p4        ->at(i));
 
	if(usingPATJets) {
	  temp_jets_pat_genParton_id      .push_back(jets_pat_genParton_id_h       ->at(i));
	  temp_jets_pat_genPartonMother_id.push_back(jets_pat_genPartonMother_id_h ->at(i));
	  temp_jets_pat_partonFlavour     .push_back(jets_pat_partonFlavour_h      ->at(i));
	  temp_jets_pat_genParton_p4      .push_back(jets_pat_genParton_p4_h       ->at(i));
	  temp_jets_pat_genPartonMother_p4.push_back(jets_pat_genPartonMother_p4_h ->at(i));
	  temp_jets_pat_genJet_p4         .push_back(jets_pat_genJet_p4_h          ->at(i));
	  temp_jets_pat_jet_p4            .push_back(jets_pat_jet_p4_h             ->at(i));
	  temp_jets_pat_noCorrF           .push_back(jets_pat_noCorrF_h            ->at(i));                  
	  temp_jets_pat_udsCorrF          .push_back(jets_pat_udsCorrF_h           ->at(i));                 
	  temp_jets_pat_gluCorrF          .push_back(jets_pat_gluCorrF_h           ->at(i));                 
	  temp_jets_pat_cCorrF            .push_back(jets_pat_cCorrF_h             ->at(i));                   
	  temp_jets_pat_bCorrF            .push_back(jets_pat_bCorrF_h             ->at(i));                   
	  temp_jets_pat_jetCharge        .push_back(jets_pat_jetCharge_h          ->at(i));               
	}

      } else {
	temp_other_jets_mc_id            .push_back(jets_mc_id           ->at(i));
	temp_other_jets_emFrac           .push_back(jets_emFrac          ->at(i));
	temp_other_jets_chFrac           .push_back(jets_chFrac          ->at(i));
	  
	temp_other_jets_mc_emEnergy      .push_back(jets_mc_emEnergy     ->at(i));
	temp_other_jets_mc_hadEnergy     .push_back(jets_mc_hadEnergy    ->at(i));
	temp_other_jets_mc_invEnergy     .push_back(jets_mc_invEnergy    ->at(i));
	temp_other_jets_mc_otherEnergy   .push_back(jets_mc_otherEnergy  ->at(i));  
	temp_other_jets_cor              .push_back(jets_cor             ->at(i));
	temp_other_jets_EMFcor           .push_back(jets_EMFcor          ->at(i));
	  
	temp_other_jets_p4               .push_back(jets_p4              ->at(i));
	temp_other_jets_mc_p4            .push_back(jets_mc_p4           ->at(i));
	temp_other_jets_mc_gp_p4         .push_back(jets_mc_gp_p4        ->at(i));

	if(usingPATJets) {
	  temp_other_jets_pat_genParton_id      .push_back(jets_pat_genParton_id_h       ->at(i));
	  temp_other_jets_pat_genPartonMother_id.push_back(jets_pat_genPartonMother_id_h ->at(i));
	  temp_other_jets_pat_partonFlavour     .push_back(jets_pat_partonFlavour_h      ->at(i));
	  temp_other_jets_pat_genParton_p4      .push_back(jets_pat_genParton_p4_h       ->at(i));
	  temp_other_jets_pat_genPartonMother_p4.push_back(jets_pat_genPartonMother_p4_h ->at(i));
	  temp_other_jets_pat_genJet_p4         .push_back(jets_pat_genJet_p4_h          ->at(i));
	  temp_other_jets_pat_jet_p4            .push_back(jets_pat_jet_p4_h             ->at(i));
	  temp_other_jets_pat_noCorrF           .push_back(jets_pat_noCorrF_h            ->at(i));                  
	  temp_other_jets_pat_udsCorrF          .push_back(jets_pat_udsCorrF_h           ->at(i));                 
	  temp_other_jets_pat_gluCorrF          .push_back(jets_pat_gluCorrF_h           ->at(i));                 
	  temp_other_jets_pat_cCorrF            .push_back(jets_pat_cCorrF_h             ->at(i));                   
	  temp_other_jets_pat_bCorrF            .push_back(jets_pat_bCorrF_h             ->at(i));                   
	  temp_other_jets_pat_jetCharge        .push_back(jets_pat_jetCharge_h          ->at(i));               
	}
	  
      }
    }
       

    //push these into the hyp_jets and hyp_other_jets vars
    hyp_jets_mc_id                 ->push_back(temp_jets_mc_id              );
    hyp_jets_emFrac                ->push_back(temp_jets_emFrac             );
    hyp_jets_chFrac                ->push_back(temp_jets_chFrac             );
    hyp_jets_mc_emEnergy           ->push_back(temp_jets_mc_emEnergy        );
    hyp_jets_mc_hadEnergy          ->push_back(temp_jets_mc_hadEnergy       );
    hyp_jets_mc_invEnergy          ->push_back(temp_jets_mc_invEnergy       );
    hyp_jets_mc_otherEnergy        ->push_back(temp_jets_mc_otherEnergy     );
    hyp_jets_cor                   ->push_back(temp_jets_cor                );
    hyp_jets_EMFcor                ->push_back(temp_jets_EMFcor             );
      
    hyp_jets_p4                    ->push_back(temp_jets_p4                  );
    hyp_jets_mc_p4                 ->push_back(temp_jets_mc_p4               );
    hyp_jets_mc_gp_p4              ->push_back(temp_jets_mc_gp_p4            );
      
    hyp_other_jets_mc_id           ->push_back(temp_other_jets_mc_id         );
    hyp_other_jets_emFrac          ->push_back(temp_other_jets_emFrac        );
    hyp_other_jets_chFrac          ->push_back(temp_other_jets_chFrac        );
    hyp_other_jets_mc_emEnergy     ->push_back(temp_other_jets_mc_emEnergy   );
    hyp_other_jets_mc_hadEnergy    ->push_back(temp_other_jets_mc_hadEnergy  );
    hyp_other_jets_mc_invEnergy    ->push_back(temp_other_jets_mc_invEnergy  );
    hyp_other_jets_mc_otherEnergy  ->push_back(temp_other_jets_mc_otherEnergy);
    hyp_other_jets_cor             ->push_back(temp_other_jets_cor           );
    hyp_other_jets_EMFcor          ->push_back(temp_other_jets_EMFcor        );
      
    hyp_other_jets_p4              ->push_back(temp_other_jets_p4            );
    hyp_other_jets_mc_p4           ->push_back(temp_other_jets_mc_p4         );
    hyp_other_jets_mc_gp_p4        ->push_back(temp_other_jets_mc_gp_p4      );

   
    if(usingPATJets) {
      hyp_jets_pat_genParton_id            ->push_back(temp_jets_pat_genParton_id             );
      hyp_jets_pat_genPartonMother_id      ->push_back(temp_jets_pat_genPartonMother_id       );
      hyp_jets_pat_partonFlavour           ->push_back(temp_jets_pat_partonFlavour            );
      hyp_jets_pat_genParton_p4            ->push_back(temp_jets_pat_genParton_p4             );
      hyp_jets_pat_genPartonMother_p4      ->push_back(temp_jets_pat_genPartonMother_p4       );
      hyp_jets_pat_genJet_p4               ->push_back(temp_jets_pat_genJet_p4                );
      hyp_jets_pat_jet_p4                  ->push_back(temp_jets_pat_jet_p4                   );
      hyp_jets_pat_noCorrF                 ->push_back(temp_jets_pat_noCorrF                );                  
      hyp_jets_pat_udsCorrF                ->push_back(temp_jets_pat_udsCorrF               );                 
      hyp_jets_pat_gluCorrF                ->push_back(temp_jets_pat_gluCorrF               );                 
      hyp_jets_pat_cCorrF                  ->push_back(temp_jets_pat_cCorrF                 );                   
      hyp_jets_pat_bCorrF                  ->push_back(temp_jets_pat_bCorrF                 );                   
      hyp_jets_pat_jetCharge              ->push_back(temp_jets_pat_jetCharge                );               
	
      hyp_other_jets_pat_genParton_id      ->push_back(temp_other_jets_pat_genParton_id       );
      hyp_other_jets_pat_genPartonMother_id->push_back(temp_other_jets_pat_genPartonMother_id );
      hyp_other_jets_pat_partonFlavour     ->push_back(temp_other_jets_pat_partonFlavour      );
      hyp_other_jets_pat_genParton_p4      ->push_back(temp_other_jets_pat_genParton_p4       );
      hyp_other_jets_pat_genPartonMother_p4->push_back(temp_other_jets_pat_genPartonMother_p4 );
      hyp_other_jets_pat_genJet_p4         ->push_back(temp_other_jets_pat_genJet_p4          );
      hyp_other_jets_pat_jet_p4            ->push_back(temp_other_jets_pat_jet_p4             );
      hyp_other_jets_pat_noCorrF           ->push_back(temp_other_jets_pat_noCorrF          );                  
      hyp_other_jets_pat_udsCorrF          ->push_back(temp_other_jets_pat_udsCorrF         );                 
      hyp_other_jets_pat_gluCorrF          ->push_back(temp_other_jets_pat_gluCorrF         );                 
      hyp_other_jets_pat_cCorrF            ->push_back(temp_other_jets_pat_cCorrF           );                   
      hyp_other_jets_pat_bCorrF            ->push_back(temp_other_jets_pat_bCorrF           );                   
      hyp_other_jets_pat_jetCharge        ->push_back(temp_other_jets_pat_jetCharge          );               
    }
     
    double hypmet = *evt_met;
    double hypmetPhi = *evt_metphi;
    METUtilities::correctMETmuons_crossedE(muon_pair,
					   hypmet, hypmetPhi, mus_e_em->at(mus_index), 
					   mus_e_had->at(mus_index),  mus_e_ho->at(mus_index) );

    //now use the expected MIP deposit
    double hypmet_MIP    = *evt_met;
    double hypmetPhi_MIP = *evt_metphi;
    METUtilities::correctMETmuons_expMIP(muon_pair,
					 hypmet_MIP, hypmetPhi_MIP );
      
    //now the s9 correction
    double hypmet_S9    = *evt_met;
    double hypmetPhi_S9 = *evt_metphi;
    METUtilities::correctMETmuons_S9E(muon_pair,
				      hypmet_S9, hypmetPhi_S9, mus_e_emS9->at(mus_index), 
				      mus_e_hadS9->at(mus_index),  mus_e_hoS9->at(mus_index) );
      
    //no calo or MIP correction
    double hypmet_nocalo = *evt_met;
    double hypmetPhi_nocalo = *evt_metphi;
    METUtilities::correctMETmuons_nocalo(muon_pair,
					 hypmet_nocalo, hypmetPhi_nocalo);
            
    double metJes5 = hypmet;
    double metPhiJes5 = hypmetPhi;
    double metJes10 = hypmet;
    double metPhiJes10 = hypmetPhi;
    double metJes15 = hypmet;
    double metPhiJes15 = hypmetPhi;
    double metJes30 = hypmet;
    double metPhiJes30 = hypmetPhi;
    double metJes50 = hypmet;
    double metPhiJes50 = hypmetPhi;
      
    //       double metEMF5 = hypmet;
    //       double metPhiEMF5 = hypmetPhi;
    //       double metEMF10 = hypmet;
    //       double metPhiEMF10 = hypmetPhi;
    //       double metEMF15 = hypmet;
    //       double metPhiEMF15 = hypmetPhi;
    //       double metEMF30 = hypmet;
    //       double metPhiEMF30 = hypmetPhi;
    //       double metEMF50 = hypmet;
    //       double metPhiEMF50 = hypmetPhi;
       
        
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes5, metPhiJes5, 5);
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes10, metPhiJes10, 10);
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes15, metPhiJes15, 15);
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes30, metPhiJes30, 30);
    METUtilities::correctedJetMET(jets_noel_p4, jets_noel_jescor,
				  metJes50, metPhiJes50, 50);
      
             
    //   METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
    // 				    metEMF5, metPhiEMF5, 5);
    //       METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
    // 				    metEMF10, metPhiEMF10, 10);
    //       METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
    // 				    metEMF15, metPhiEMF15, 15);
    //       METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
    // 				    metEMF30, metPhiEMF30, 30);
    //       METUtilities::correctedJetMET(jets_noel_p4, jets_noel_EMFcor,
    // 				    metEMF50, metPhiEMF50, 50);

      
    hyp_njets         ->push_back(temp_jets_p4.size()                   );     
    hyp_nojets        ->push_back(temp_other_jets_p4.size()             );     
    hyp_p4            ->push_back(mus_p4->at(mus_index)+
				  els_p4->at(els_index)                 );
	
	
	
    if(el_pt < tightptcut && mu_pt > tightptcut) {
      hyp_type            ->push_back(1);
	  
      hyp_lt_validHits    ->push_back(mus_validHits    ->at(mus_index)  );
      hyp_lt_lostHits     ->push_back(mus_lostHits     ->at(mus_index)  );
      hyp_lt_mc_id        ->push_back(mus_mc_id        ->at(mus_index)  );
      hyp_lt_charge       ->push_back(mus_charge       ->at(mus_index)  );
      hyp_lt_mc_motherid  ->push_back(mus_mc_motherid  ->at(mus_index)  );
      hyp_lt_index        ->push_back(mus_index                         );
      hyp_lt_id           ->push_back(-13*(mus_charge   ->at(mus_index)) );
      hyp_lt_d0           ->push_back(mus_d0           ->at(mus_index)  );
      hyp_lt_z0           ->push_back(mus_z0           ->at(mus_index)  );
      hyp_lt_d0corr       ->push_back(mus_d0corr       ->at(mus_index)  );
      hyp_lt_z0corr       ->push_back(mus_z0corr       ->at(mus_index)  );
      hyp_lt_vertexphi    ->push_back(mus_vertexphi    ->at(mus_index)  );
      hyp_lt_chi2         ->push_back(mus_chi2         ->at(mus_index)  );
      hyp_lt_ndof         ->push_back(mus_ndof         ->at(mus_index)  );
      hyp_lt_d0Err        ->push_back(mus_d0Err        ->at(mus_index)  );
      hyp_lt_z0Err        ->push_back(mus_z0Err        ->at(mus_index)  );
      hyp_lt_ptErr        ->push_back(mus_ptErr        ->at(mus_index)  );
      hyp_lt_etaErr       ->push_back(mus_etaErr       ->at(mus_index)  );
      hyp_lt_phiErr       ->push_back(mus_phiErr       ->at(mus_index)  );
      hyp_lt_outerPhi     ->push_back(mus_outerPhi     ->at(mus_index)  );
      hyp_lt_outerEta     ->push_back(mus_outerEta     ->at(mus_index)  );
      hyp_lt_iso          ->push_back(mus_iso          ->at(mus_index)  );
      hyp_lt_tkIso        ->push_back(mus_iso03_sumPt  ->at(mus_index)  );
      hyp_lt_p4           ->push_back(mus_p4           ->at(mus_index)  );
      hyp_lt_trk_p4       ->push_back(mus_trk_p4       ->at(mus_index)  );
      hyp_lt_mc_p4        ->push_back(mus_mc_p4        ->at(mus_index)  );

      hyp_ll_validHits    ->push_back(els_validHits    ->at(els_index)  );
      hyp_ll_lostHits     ->push_back(els_lostHits     ->at(els_index)  );
      hyp_ll_mc_id        ->push_back(els_mc_id        ->at(els_index)  );
      hyp_ll_charge       ->push_back(els_charge       ->at(els_index)  );
      hyp_ll_mc_motherid  ->push_back(els_mc_motherid  ->at(els_index)  );
      hyp_ll_index        ->push_back(els_index                         );
      hyp_ll_id           ->push_back(-11*(els_charge   ->at(els_index)) );
      hyp_ll_d0           ->push_back(els_d0           ->at(els_index)  );
      hyp_ll_z0           ->push_back(els_z0           ->at(els_index)  );
      hyp_ll_d0corr       ->push_back(els_d0corr       ->at(els_index)  );
      hyp_ll_z0corr       ->push_back(els_z0corr       ->at(els_index)  );
      hyp_ll_vertexphi    ->push_back(els_vertexphi    ->at(els_index)  );
      hyp_ll_chi2         ->push_back(els_chi2         ->at(els_index)  );
      hyp_ll_ndof         ->push_back(els_ndof         ->at(els_index)  );
      hyp_ll_d0Err        ->push_back(els_d0Err        ->at(els_index)  );
      hyp_ll_z0Err        ->push_back(els_z0Err        ->at(els_index)  );
      hyp_ll_ptErr        ->push_back(els_ptErr        ->at(els_index)  );
      hyp_ll_etaErr       ->push_back(els_etaErr       ->at(els_index)  );
      hyp_ll_phiErr       ->push_back(els_phiErr       ->at(els_index)  );
      hyp_ll_outerPhi     ->push_back(els_outerPhi     ->at(els_index)  );
      hyp_ll_outerEta     ->push_back(els_outerEta     ->at(els_index)  );
      hyp_ll_iso          ->push_back(els_tkIso        ->at(els_index)  );
      hyp_ll_tkIso        ->push_back(els_tkIso        ->at(els_index)  );
      hyp_ll_p4           ->push_back(els_p4           ->at(els_index)  );
      hyp_ll_trk_p4       ->push_back(els_trk_p4       ->at(els_index)  );
      hyp_ll_mc_p4        ->push_back(els_mc_p4        ->at(els_index)  );
	  
	  
    } else {
      hyp_type            ->push_back(2);
		  
      hyp_lt_validHits    ->push_back(els_validHits    ->at(els_index)  );
      hyp_lt_lostHits     ->push_back(els_lostHits     ->at(els_index)  );
      hyp_lt_mc_id        ->push_back(els_mc_id        ->at(els_index)  );
      hyp_lt_charge       ->push_back(els_charge       ->at(els_index)  );
      hyp_lt_mc_motherid  ->push_back(els_mc_motherid  ->at(els_index)  );
      hyp_lt_index        ->push_back(els_index                         );
      hyp_lt_id           ->push_back(-11*(els_charge   ->at(els_index)) );
      hyp_lt_d0           ->push_back(els_d0           ->at(els_index)  );
      hyp_lt_z0           ->push_back(els_z0           ->at(els_index)  );
      hyp_lt_d0corr       ->push_back(els_d0corr       ->at(els_index)  );
      hyp_lt_z0corr       ->push_back(els_z0corr       ->at(els_index)  );
      hyp_lt_vertexphi    ->push_back(els_vertexphi    ->at(els_index)  );
      hyp_lt_chi2         ->push_back(els_chi2         ->at(els_index)  );
      hyp_lt_ndof         ->push_back(els_ndof         ->at(els_index)  );
      hyp_lt_d0Err        ->push_back(els_d0Err        ->at(els_index)  );
      hyp_lt_z0Err        ->push_back(els_z0Err        ->at(els_index)  );
      hyp_lt_ptErr        ->push_back(els_ptErr        ->at(els_index)  );
      hyp_lt_etaErr       ->push_back(els_etaErr       ->at(els_index)  );
      hyp_lt_phiErr       ->push_back(els_phiErr       ->at(els_index)  );
      hyp_lt_outerPhi     ->push_back(els_outerPhi     ->at(els_index)  );
      hyp_lt_outerEta     ->push_back(els_outerEta     ->at(els_index)  );
      hyp_lt_iso          ->push_back(els_tkIso        ->at(els_index)  );
      hyp_lt_tkIso        ->push_back(els_tkIso        ->at(els_index)  );
      hyp_lt_p4           ->push_back(els_p4           ->at(els_index)  );
      hyp_lt_trk_p4       ->push_back(els_trk_p4       ->at(els_index)  );
      hyp_lt_mc_p4        ->push_back(els_mc_p4        ->at(els_index)  );
	   

      
      hyp_ll_validHits    ->push_back(mus_validHits    ->at(mus_index)  );
      hyp_ll_lostHits     ->push_back(mus_lostHits     ->at(mus_index)  );
      hyp_ll_mc_id        ->push_back(mus_mc_id        ->at(mus_index)  );
      hyp_ll_charge       ->push_back(mus_charge       ->at(mus_index)  );
      hyp_ll_mc_motherid  ->push_back(mus_mc_motherid  ->at(mus_index)  );
      hyp_ll_index        ->push_back(mus_index                         );
      hyp_ll_id           ->push_back(-13*(mus_charge   ->at(mus_index)) );
      hyp_ll_d0           ->push_back(mus_d0           ->at(mus_index)  );
      hyp_ll_z0           ->push_back(mus_z0           ->at(mus_index)  );
      hyp_ll_d0corr       ->push_back(mus_d0corr       ->at(mus_index)  );
      hyp_ll_z0corr       ->push_back(mus_z0corr       ->at(mus_index)  );
      hyp_ll_vertexphi    ->push_back(mus_vertexphi    ->at(mus_index)  );
      hyp_ll_chi2         ->push_back(mus_chi2         ->at(mus_index)  );
      hyp_ll_ndof         ->push_back(mus_ndof         ->at(mus_index)  );
      hyp_ll_d0Err        ->push_back(mus_d0Err        ->at(mus_index)  );
      hyp_ll_z0Err        ->push_back(mus_z0Err        ->at(mus_index)  );
      hyp_ll_ptErr        ->push_back(mus_ptErr        ->at(mus_index)  );
      hyp_ll_etaErr       ->push_back(mus_etaErr       ->at(mus_index)  );
      hyp_ll_phiErr       ->push_back(mus_phiErr       ->at(mus_index)  );
      hyp_ll_outerPhi     ->push_back(mus_outerPhi     ->at(mus_index)  );
      hyp_ll_outerEta     ->push_back(mus_outerEta     ->at(mus_index)  );
      hyp_ll_iso          ->push_back(mus_iso          ->at(mus_index)  );
      hyp_ll_tkIso        ->push_back(mus_iso03_sumPt  ->at(mus_index)  );
      hyp_ll_p4           ->push_back(mus_p4           ->at(mus_index)  );
      hyp_ll_trk_p4       ->push_back(mus_trk_p4       ->at(mus_index)  );
      hyp_ll_mc_p4        ->push_back(mus_mc_p4        ->at(mus_index)  );

    }

    hyp_met             ->push_back(hypmet                                 );
    hyp_metPhi          ->push_back(hypmetPhi                              );
    hyp_metCaloExp      ->push_back(hypmet_MIP                             );
    hyp_metPhiCaloExp   ->push_back(hypmetPhi_MIP                          );
    hyp_metCone         ->push_back(hypmet_S9                              );
    hyp_metPhiCone      ->push_back(hypmetPhi_S9                           );
    hyp_metNoCalo       ->push_back(hypmet_nocalo                          );
    hyp_metPhiNoCalo    ->push_back(hypmetPhi_nocalo                       );
    hyp_metAll          ->push_back(metAll                                 );
    hyp_metPhiAll       ->push_back(metPhiAll                              );
    hyp_metAllCaloExp   ->push_back(metAllCaloExp                          );
    hyp_metPhiAllCaloExp->push_back(metPhiAllCaloExp                       );
    hyp_metJes5         ->push_back(metJes5                                );
    hyp_metPhiJes5      ->push_back(metPhiJes5                             );
    hyp_metJes10        ->push_back(metJes10                               );
    hyp_metPhiJes10     ->push_back(metPhiJes10                            );
    hyp_metJes15        ->push_back(metJes15                               );
    hyp_metPhiJes15     ->push_back(metPhiJes15                            );
    hyp_metJes30        ->push_back(metJes30                               );
    hyp_metPhiJes30     ->push_back(metPhiJes30                            );
    hyp_metJes50        ->push_back(metJes50                               );
    hyp_metPhiJes50     ->push_back(metPhiJes50                            );
    // hyp_metEMF5         ->push_back(metEMF5                                );
    //       hyp_metPhiEMF5      ->push_back(metPhiEMF5                             );
    //       hyp_metEMF10        ->push_back(metEMF10                               );
    //       hyp_metPhiEMF10     ->push_back(metPhiEMF10                            );
    //       hyp_metEMF15        ->push_back(metEMF15                               );
    //       hyp_metPhiEMF15     ->push_back(metPhiEMF15                            );
    //       hyp_metEMF30        ->push_back(metEMF30                               );
    //       hyp_metPhiEMF30     ->push_back(metPhiEMF30                            );
    //       hyp_metEMF50        ->push_back(metEMF50                               );
    //       hyp_metPhiEMF50     ->push_back(metPhiEMF50                            );
    hyp_metDPhiJet10    ->push_back(METUtilities::metObjDPhi(*jets_p4,
							     hypmetPhi,
							     10)           );
    hyp_metDPhiJet15    ->push_back(METUtilities::metObjDPhi(*jets_p4,
							     hypmetPhi,
							     15)           ); 
    hyp_metDPhiJet20    ->push_back(METUtilities::metObjDPhi(*jets_p4,
							     hypmetPhi,
							     20)           );
    hyp_metDPhiTrk10    ->push_back(METUtilities::metObjDPhi(*trks_p4,
							     hypmetPhi,
							     10)           );
    hyp_metDPhiTrk25    ->push_back(METUtilities::metObjDPhi(*trks_p4,
							     hypmetPhi,
							     25)           );
    hyp_metDPhiTrk50    ->push_back(METUtilities::metObjDPhi(*trks_p4,
							     hypmetPhi,
							     50)           );

		
  }
 }
  



  
 iEvent.put(hyp_type                     ,"hyptype"                     );
 iEvent.put(hyp_njets                    ,"hypnjets"                    );
 iEvent.put(hyp_nojets                   ,"hypnojets"                   );
 iEvent.put(hyp_p4                      ,"hypp4"                        );
 
 iEvent.put(hyp_lt_validHits             ,"hypltvalidHits"              );
 iEvent.put(hyp_lt_lostHits              ,"hypltlostHits"               );
 iEvent.put(hyp_lt_mc_id                 ,"hypltmcid"                   );
 iEvent.put(hyp_lt_charge                ,"hypltcharge"                 );
 iEvent.put(hyp_lt_mc_motherid           ,"hypltmcmotherid"             );
 iEvent.put(hyp_lt_index                 ,"hypltindex"                  );
 iEvent.put(hyp_lt_id                    ,"hypltid"                     );
 iEvent.put(hyp_lt_d0                    ,"hypltd0"                     );
 iEvent.put(hyp_lt_z0                    ,"hypltz0"                     );
 iEvent.put(hyp_lt_d0corr                ,"hypltd0corr"                 );
 iEvent.put(hyp_lt_z0corr                ,"hypltz0corr"                 );
 iEvent.put(hyp_lt_vertexphi             ,"hypltvertexphi"              );
 iEvent.put(hyp_lt_chi2                  ,"hypltchi2"                   );
 iEvent.put(hyp_lt_ndof                  ,"hypltndof"                   );
 iEvent.put(hyp_lt_d0Err                 ,"hypltd0Err"                  );
 iEvent.put(hyp_lt_z0Err                 ,"hypltz0Err"                  );
 iEvent.put(hyp_lt_ptErr                 ,"hypltptErr"                  );
 iEvent.put(hyp_lt_etaErr                ,"hypltetaErr"                 );
 iEvent.put(hyp_lt_phiErr                ,"hypltphiErr"                 );
 iEvent.put(hyp_lt_outerPhi              ,"hypltouterPhi"               );
 iEvent.put(hyp_lt_outerEta              ,"hypltouterEta"               );
 iEvent.put(hyp_lt_iso                   ,"hypltiso"                    );
 iEvent.put(hyp_lt_tkIso                 ,"hyplttkIso"                  );
 iEvent.put(hyp_lt_p4                    ,"hypltp4"                     );
 iEvent.put(hyp_lt_trk_p4                ,"hyplttrkp4"                  );
 iEvent.put(hyp_lt_mc_p4                 ,"hypltmcp4"                   );
 
 iEvent.put(hyp_ll_validHits             ,"hypllvalidHits"              );
 iEvent.put(hyp_ll_lostHits              ,"hyplllostHits"               );
 iEvent.put(hyp_ll_mc_id                 ,"hypllmcid"                   );
 iEvent.put(hyp_ll_charge                ,"hypllcharge"                 );
 iEvent.put(hyp_ll_mc_motherid           ,"hypllmcmotherid"             );
 iEvent.put(hyp_ll_index                 ,"hypllindex"                  );
 iEvent.put(hyp_ll_id                    ,"hypllid"                     );
 iEvent.put(hyp_ll_d0                    ,"hyplld0"                     );
 iEvent.put(hyp_ll_z0                    ,"hypllz0"                     );
 iEvent.put(hyp_ll_d0corr                ,"hyplld0corr"                 );
 iEvent.put(hyp_ll_z0corr                ,"hypllz0corr"                 );
 iEvent.put(hyp_ll_vertexphi             ,"hypllvertexphi"              );
 iEvent.put(hyp_ll_chi2                  ,"hypllchi2"                   );
 iEvent.put(hyp_ll_ndof                  ,"hypllndof"                   );
 iEvent.put(hyp_ll_d0Err                 ,"hyplld0Err"                  );
 iEvent.put(hyp_ll_z0Err                 ,"hypllz0Err"                  );
 iEvent.put(hyp_ll_ptErr                 ,"hypllptErr"                  );
 iEvent.put(hyp_ll_etaErr                ,"hyplletaErr"                 );
 iEvent.put(hyp_ll_phiErr                ,"hypllphiErr"                 );
 iEvent.put(hyp_ll_outerPhi              ,"hypllouterPhi"               );
 iEvent.put(hyp_ll_outerEta              ,"hypllouterEta"               );
 iEvent.put(hyp_ll_iso                   ,"hyplliso"                    );
 iEvent.put(hyp_ll_tkIso                 ,"hyplltkIso"                  );
 iEvent.put(hyp_ll_p4                    ,"hypllp4"                     );
 iEvent.put(hyp_ll_trk_p4                ,"hyplltrkp4"                  );
 iEvent.put(hyp_ll_mc_p4                 ,"hypllmcp4"                   );
 
 iEvent.put(hyp_met                      ,"hypmet"                      );
 iEvent.put(hyp_metPhi                   ,"hypmetPhi"                   );
 iEvent.put(hyp_metCaloExp               ,"hypmetCaloExp"               );
 iEvent.put(hyp_metPhiCaloExp            ,"hypmetPhiCaloExp"            );
 iEvent.put(hyp_metCone                  ,"hypmetCone"                  );
 iEvent.put(hyp_metPhiCone               ,"hypmetPhiCone"               );
 iEvent.put(hyp_metNoCalo                ,"hypmetNoCalo"                );
 iEvent.put(hyp_metPhiNoCalo             ,"hypmetPhiNoCalo"             );
 iEvent.put(hyp_metAll                   ,"hypmetAll"                   );
 iEvent.put(hyp_metPhiAll                ,"hypmetPhiAll"                );
 iEvent.put(hyp_metAllCaloExp            ,"hypmetAllCaloExp"            );
 iEvent.put(hyp_metPhiAllCaloExp         ,"hypmetPhiAllCaloExp"         );
 iEvent.put(hyp_metJes5                  ,"hypmetJes5"                  );
 iEvent.put(hyp_metPhiJes5               ,"hypmetPhiJes5"               );
 iEvent.put(hyp_metJes10                 ,"hypmetJes10"                 );
 iEvent.put(hyp_metPhiJes10              ,"hypmetPhiJes10"              );
 iEvent.put(hyp_metJes15                 ,"hypmetJes15"                 );
 iEvent.put(hyp_metPhiJes15              ,"hypmetPhiJes15"              );
 iEvent.put(hyp_metJes30                 ,"hypmetJes30"                 );
 iEvent.put(hyp_metPhiJes30              ,"hypmetPhiJes30"              );
 iEvent.put(hyp_metJes50                 ,"hypmetJes50"                 );
 iEvent.put(hyp_metPhiJes50              ,"hypmetPhiJes50"              );
 //   iEvent.put(hyp_metEMF5                  ,"hypmetEMF5"                  );
//   iEvent.put(hyp_metPhiEMF5               ,"hypmetPhiEMF5"               );
//   iEvent.put(hyp_metEMF10                 ,"hypmetEMF10"                 );
//   iEvent.put(hyp_metPhiEMF10              ,"hypmetPhiEMF10"              );
//   iEvent.put(hyp_metEMF15                 ,"hypmetEMF15"                 );
//   iEvent.put(hyp_metPhiEMF15              ,"hypmetPhiEMF15"              );
//   iEvent.put(hyp_metEMF30                 ,"hypmetEMF30"                 );
//   iEvent.put(hyp_metPhiEMF30              ,"hypmetPhiEMF30"              );
//   iEvent.put(hyp_metEMF50                 ,"hypmetEMF50"                 );
//   iEvent.put(hyp_metPhiEMF50              ,"hypmetPhiEMF50"              );
iEvent.put(hyp_metDPhiJet10             ,"hypmetDPhiJet10"             );
iEvent.put(hyp_metDPhiJet15             ,"hypmetDPhiJet15"             );
iEvent.put(hyp_metDPhiJet20             ,"hypmetDPhiJet20"             );
iEvent.put(hyp_metDPhiTrk10             ,"hypmetDPhiTrk10"             );
iEvent.put(hyp_metDPhiTrk25             ,"hypmetDPhiTrk25"             );
iEvent.put(hyp_metDPhiTrk50             ,"hypmetDPhiTrk50"             );
  
  
iEvent.put(hyp_jets_mc_id               ,"hypjetsmcid"                 );
iEvent.put(hyp_jets_emFrac              ,"hypjetsemFrac"               );
iEvent.put(hyp_jets_chFrac              ,"hypjetschFrac"               );
iEvent.put(hyp_jets_mc_emEnergy         ,"hypjetsmcemEnergy"           );
iEvent.put(hyp_jets_mc_hadEnergy        ,"hypjetsmchadEnergy"          );
iEvent.put(hyp_jets_mc_invEnergy        ,"hypjetsmcinvEnergy"          );
iEvent.put(hyp_jets_mc_otherEnergy      ,"hypjetsmcotherEnergy"        );
iEvent.put(hyp_jets_cor                 ,"hypjetscor"                  );
iEvent.put(hyp_jets_EMFcor              ,"hypjetsEMFcor"               );
iEvent.put(hyp_other_jets_mc_id         ,"hypotherjetsmcid"            );
iEvent.put(hyp_other_jets_emFrac        ,"hypotherjetsemFrac"          );
iEvent.put(hyp_other_jets_chFrac        ,"hypotherjetschFrac"          );
iEvent.put(hyp_other_jets_mc_emEnergy   ,"hypotherjetsmcemEnergy"      );
iEvent.put(hyp_other_jets_mc_hadEnergy  ,"hypotherjetsmchadEnergy"     );
iEvent.put(hyp_other_jets_mc_invEnergy  ,"hypotherjetsmcinvEnergy"     );
iEvent.put(hyp_other_jets_mc_otherEnergy,"hypotherjetsmcotherEnergy"   );
iEvent.put(hyp_other_jets_cor           ,"hypotherjetscor"             );
iEvent.put(hyp_other_jets_EMFcor        ,"hypotherjetsEMFcor"          );
iEvent.put(hyp_jets_p4                  ,"hypjetsp4"                   );         
iEvent.put(hyp_jets_mc_p4               ,"hypjetsmcp4"                 );      
iEvent.put(hyp_jets_mc_gp_p4            ,"hypjetsmcgpp4"               );   
iEvent.put(hyp_other_jets_p4            ,"hypotherjetsp4"              );      
iEvent.put(hyp_other_jets_mc_p4         ,"hypotherjetsmcp4"            );   
iEvent.put(hyp_other_jets_mc_gp_p4      ,"hypotherjetsmcgpp4"          );


if(usingPATJets) {
  iEvent.put(hyp_jets_pat_genParton_id             ,"hypjetspatgenPartonid"              );    //ok            
  iEvent.put(hyp_jets_pat_genPartonMother_id       ,"hypjetspatgenPartonMotherid"        );    //ok      
  iEvent.put(hyp_jets_pat_partonFlavour            ,"hypjetspatpartonFlavour"             );   //ok            
  iEvent.put(hyp_other_jets_pat_genParton_id       ,"hypotherjetspatgenPartonid"        );     //ok      
  iEvent.put(hyp_other_jets_pat_genPartonMother_id ,"hypotherjetspatgenPartonMotherid"  );    
  iEvent.put(hyp_other_jets_pat_partonFlavour      ,"hypotherjetspatpartonFlavour"       );         
  
  iEvent.put(hyp_jets_pat_genParton_p4             ,"hypjetspatgenPartonp4"              );      
  iEvent.put(hyp_jets_pat_genPartonMother_p4       ,"hypjetspatgenPartonMotherp4"        );
  iEvent.put(hyp_jets_pat_genJet_p4                ,"hypjetspatgenJetp4"                 );         
  iEvent.put(hyp_jets_pat_jet_p4                   ,"hypjetspatjetp4"                    );            
  iEvent.put(hyp_other_jets_pat_genParton_p4       ,"hypotherjetspatgenPartonp4"        );
  iEvent.put(hyp_other_jets_pat_genPartonMother_p4 ,"hypotherjetspatgenPartonMotherp4"  );
  iEvent.put(hyp_other_jets_pat_genJet_p4          ,"hypotherjetspatgenJetp4"           );   
  iEvent.put(hyp_other_jets_pat_jet_p4             ,"hypotherjetspatjetp4"              );      
 }  

  
}

// ------------ method called once each job just before starting event loop  ------------
void HypDilepMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HypDilepMaker::endJob() {
}

//----------------------------------------------------------------------------------------
bool HypDilepMaker::testJetForElectrons(const LorentzVector& jetP4, const LorentzVector& elP4) {
  
  
  bool matched = false;
  float elphi  = elP4.Phi();
  float jetphi = jetP4.Phi();
   
  float eleta  = elP4.Eta();
  float jeteta = jetP4.Eta();
   
  float dphi = elphi - jetphi;
  float deta = eleta - jeteta;
  if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
   
  double dR = sqrt(dphi*dphi + deta*deta);
  if (dR < 0.4) 
    matched = true;
  
  return !matched;
}


//define this as a plug-in
DEFINE_FWK_MODULE(HypDilepMaker);

