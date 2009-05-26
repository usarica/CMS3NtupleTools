//-*- C++ -*-
//
// Package:    PFTauMaker
// Class:      PFTauMaker
// 
/**\class PFTauMaker PFTauMaker.cc CMS2/NtupleMaker/src/PFTauMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// $Id: PFTauMaker.cc,v 1.1 2009/05/26 23:24:06 yanjuntu Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/PFTauMaker.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"


typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

PFTauMaker::PFTauMaker(const edm::ParameterSet& iConfig) {

  produces<vector<LorentzVector> >  ("tauspfp4"                            ).setBranchAlias("taus_pf_p4"                             );
  produces<vector<int> >            ("tauspfcharge"                        ).setBranchAlias("taus_pf_charge"                         );
  produces<vector<int> >            ("tauspfsignchargecand"                ).setBranchAlias("taus_pf_sig_ncharge_cand"               );
  produces<vector<int> >            ("tauspfisonchargecand"                ).setBranchAlias("taus_pf_iso_ncharge_cand"               );
  produces<vector<int> >            ("tauspfsignneutrcand"                 ).setBranchAlias("taus_pf_sig_nneutr_cand"                );
  produces<vector<int> >            ("tauspfisonneutrcand"                 ).setBranchAlias("taus_pf_iso_nneutr_cand"                );
  produces<vector<int> >            ("tauspfsigngammacand"                 ).setBranchAlias("taus_pf_sig_ngamma_cand"                );
  produces<vector<int> >            ("tauspfisongammacand"                 ).setBranchAlias("taus_pf_iso_ngamma_cand"                );
  produces<vector<LorentzVector> >  ("tauspfleadchargecandp4"              ).setBranchAlias("taus_pf_lead_chargecand_p4"             );
  produces<vector<LorentzVector> >  ("tauspfleadneutrcandp4"               ).setBranchAlias("taus_pf_lead_neutrcand_p4"              );
  produces<vector<float> >          ("tauspfleadchargecanSignedSipt"       ).setBranchAlias("taus_pf_lead_chargecand_Signed_Sipt"    );
  produces<vector<float> >          ("tauspfisolationchargecandPtSum"      ).setBranchAlias("taus_pf_isolationchargecandPtSum"       ); 
  produces<vector<float> >          ("tauspfisolationgammacandEtSum"       ).setBranchAlias("taus_pf_isolationgammacandEtSum"        ); 
  produces<vector<float> >          ("tauspfmaximumHCALPFClusterEt"        ).setBranchAlias("taus_pf_maximumHCALPFClusterEt"         ); 
  //electron rejection
  produces<vector<float> >          ("tauspfemf"                           ).setBranchAlias("taus_pf_emf"                            ); 
  produces<vector<float> >          ("tauspfhcalTotOverPLead"              ).setBranchAlias("taus_pf_hcalTotOverPLead"               ); 
  produces<vector<float> >          ("tauspfhcalMaxOverPLead"              ).setBranchAlias("taus_pf_hcalMaxOverPLead"               ); 
  produces<vector<float> >          ("tauspfhcal3x3OverPLead"              ).setBranchAlias("taus_pf_hcal3x3OverPLead"               ); 
  produces<vector<float> >          ("tauspfecalStripSumEOverPLead"        ).setBranchAlias("taus_pf_ecalStripSumEOverPLead"         ); 
  produces<vector<float> >          ("tauspfbremsRecoveryEOverPLead"       ).setBranchAlias("taus_pf_bremsRecoveryEOverPLead"        ); 
  produces<vector<int> >            ("tauspfelectronPreID"                 ).setBranchAlias("taus_pf_electronPreID"                  ); 
  //muon rejection
  produces<vector<int> >            ("tauspfmuonPreID"                     ).setBranchAlias("taus_pf_muonPreID"                      ); 
  produces<vector<int> >            ("tauspfhasMuonReference"              ).setBranchAlias("taus_pf_hasMuonReference"               ); 
  produces<vector<float> >          ("tauspfcaloComp"                      ).setBranchAlias("taus_pf_caloComp"                       ); 
  produces<vector<float> >          ("tauspfsegComp"                       ).setBranchAlias("taus_pf_segComp"                        ); 

  produces<vector<LorentzVector> >  ("tauspfleadtrkp4"                     ).setBranchAlias("taus_pf_leadtrk_p4"                     );
  produces<vector<float> >          ("tauspfleadtrkd0"                     ).setBranchAlias("taus_pf_leadtrk_d0"                     );  
  produces<vector<float> >          ("tauspfleadtrkz0"                     ).setBranchAlias("taus_pf_leadtrk_z0"                     ); 
  produces<vector<float> >          ("tauspfleadtrkchi2"                   ).setBranchAlias("taus_pf_leadtrk_chi2"                   ); 
  produces<vector<float> >          ("tauspfleadtrkndof"                   ).setBranchAlias("taus_pf_leadtrk_ndof"                   ); 
  produces<vector<float> >          ("tauspfleadtrkvalidHits"              ).setBranchAlias("taus_pf_leadtrk_validHits"              ); 
  produces<vector<float> >          ("tauspfleadtrklostHits"               ).setBranchAlias("taus_pf_leadtrk_lostHits"               ); 

 
 
  
  
    
   
//get setup parameters
  pftausInputTag      = iConfig.getParameter<InputTag>("pftausInputTag");

}


PFTauMaker::~PFTauMaker() {}

void  PFTauMaker::beginJob(const edm::EventSetup&) {
}

void PFTauMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void PFTauMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<LorentzVector> > taus_pf_p4                              (new vector<LorentzVector>) ;
  auto_ptr<vector<int> >           taus_pf_charge                          (new vector<int>) ;
  auto_ptr<vector<int> >           taus_pf_sig_ncharge_cand                (new vector<int>) ;
  auto_ptr<vector<int> >           taus_pf_iso_ncharge_cand                (new vector<int>) ;
  auto_ptr<vector<int> >           taus_pf_sig_nneutr_cand                 (new vector<int>) ;
  auto_ptr<vector<int> >           taus_pf_iso_nneutr_cand                 (new vector<int>) ;
  auto_ptr<vector<int> >           taus_pf_sig_ngamma_cand                 (new vector<int>) ;
  auto_ptr<vector<int> >           taus_pf_iso_ngamma_cand                 (new vector<int>) ;
  auto_ptr<vector<LorentzVector> > taus_pf_lead_chargecand_p4              (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > taus_pf_lead_neutrcand_p4               (new vector<LorentzVector>) ;
 
  auto_ptr<vector<float> >         taus_pf_lead_chargecand_Signed_Sipt     (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_isolationchargecandPtSum        (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_isolationgammacandEtSum         (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_maximumHCALPFClusterEt          (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_emf                             (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_hcalTotOverPLead                (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_hcalMaxOverPLead                (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_hcal3x3OverPLead                (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_ecalStripSumEOverPLead          (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_bremsRecoveryEOverPLead         (new vector<float>) ;
 
  auto_ptr<vector<int> >           taus_pf_electronPreID                   (new vector<int>) ;
  auto_ptr<vector<int> >           taus_pf_muonPreID                       (new vector<int>) ;
  auto_ptr<vector<int> >           taus_pf_hasMuonReference                (new vector<int>) ;
  auto_ptr<vector<float> >         taus_pf_caloComp                        (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_segComp                         (new vector<float>) ;
  auto_ptr<vector<LorentzVector> > taus_pf_leadtrk_p4                      (new vector<LorentzVector>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_d0                      (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_z0                      (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_chi2                    (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_ndof                    (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_validHits               (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_lostHits                (new vector<float>) ;
 

 
  
  Handle<View<reco::PFTau> > taus_pf_h;
  iEvent.getByLabel(pftausInputTag, taus_pf_h);
  
 
  
 for(View<reco::PFTau>::const_iterator tau_pf = taus_pf_h->begin();
      tau_pf != taus_pf_h->end(); tau_pf++) {
   if(tau_pf->leadPFChargedHadrCand().isNull()) continue;
   if(tau_pf->leadPFChargedHadrCand()->pt()<5.0) continue;
   //printf("%s  \n", "PFTau  ");

   taus_pf_p4                               ->push_back( tau_pf->p4()                                       );
   taus_pf_charge                           ->push_back( tau_pf->charge()                                   );
   taus_pf_sig_ncharge_cand                 ->push_back( tau_pf->signalPFChargedHadrCands().size()          ); 
   taus_pf_iso_ncharge_cand                 ->push_back( tau_pf->isolationPFChargedHadrCands().size()       );                                                                            
   taus_pf_sig_nneutr_cand                  ->push_back( tau_pf->signalPFNeutrHadrCands().size()            ); 
   taus_pf_iso_nneutr_cand                  ->push_back( tau_pf->isolationPFNeutrHadrCands().size()         ); 
   taus_pf_sig_ngamma_cand                  ->push_back( tau_pf->signalPFGammaCands().size()                ); 
   taus_pf_iso_ngamma_cand                  ->push_back( tau_pf->isolationPFGammaCands().size()             ); 
  
 
   taus_pf_lead_chargecand_p4               ->push_back( tau_pf->leadPFChargedHadrCand().get()->p4()        );
   if(tau_pf->leadPFNeutralCand().isNonnull())
   taus_pf_lead_neutrcand_p4                ->push_back( tau_pf->leadPFNeutralCand().get()->p4()            );
    taus_pf_lead_chargecand_Signed_Sipt      ->push_back( tau_pf->leadPFChargedHadrCandsignedSipt()); 
   taus_pf_isolationchargecandPtSum         ->push_back( tau_pf->isolationPFChargedHadrCandsPtSum() ); 
   taus_pf_isolationgammacandEtSum          ->push_back( tau_pf->isolationPFGammaCandsEtSum()       ); 
   taus_pf_maximumHCALPFClusterEt           ->push_back( tau_pf->maximumHCALPFClusterEt()       ); 
   taus_pf_emf                              ->push_back( tau_pf->emFraction()                       ); 
   taus_pf_hcalTotOverPLead                 ->push_back( tau_pf->hcalTotOverPLead()                 ); 
   taus_pf_hcalMaxOverPLead                 ->push_back( tau_pf->hcalMaxOverPLead()                 ); 
   taus_pf_hcal3x3OverPLead                 ->push_back( tau_pf->hcal3x3OverPLead()                 ); 
   taus_pf_ecalStripSumEOverPLead           ->push_back( tau_pf->ecalStripSumEOverPLead()           ); 
   taus_pf_bremsRecoveryEOverPLead          ->push_back( tau_pf->bremsRecoveryEOverPLead()          ); 
 
  
   if(tau_pf->electronPreIDDecision()) taus_pf_electronPreID      ->push_back(1);
   else                                taus_pf_electronPreID      ->push_back(0);
   if(!tau_pf->muonDecision())         taus_pf_muonPreID          ->push_back(1);
   else                                taus_pf_muonPreID          ->push_back(0);
  
   if(tau_pf->hasMuonReference()){
                                      taus_pf_hasMuonReference   ->push_back(1);
                                      taus_pf_caloComp           ->push_back(tau_pf->caloComp());
                                      taus_pf_segComp            ->push_back(tau_pf->segComp());
     
   }
   else {
                                      taus_pf_hasMuonReference   ->push_back(0);
                                      taus_pf_caloComp           ->push_back(-999.);
                                      taus_pf_segComp            ->push_back(-999.);
				     
   }
  
   const TrackRef leadTrack = tau_pf->leadTrack()  ;
   taus_pf_leadtrk_p4                    ->push_back( LorentzVector( leadTrack.get()->px(), leadTrack.get()->py(),
						leadTrack.get()->pz(), leadTrack.get()->p())                );
   taus_pf_leadtrk_d0                    ->push_back( leadTrack->d0()                                       );
   taus_pf_leadtrk_z0                    ->push_back( leadTrack->dz()                                       );
   taus_pf_leadtrk_chi2                  ->push_back( leadTrack->chi2()                                     );
   taus_pf_leadtrk_ndof                  ->push_back( leadTrack->ndof()                                     );
   taus_pf_leadtrk_validHits             ->push_back( leadTrack->numberOfValidHits()                        );
   taus_pf_leadtrk_lostHits              ->push_back( leadTrack->numberOfLostHits()                         );
  
   
 }


 iEvent.put(taus_pf_p4                                   ,"tauspfp4"                                       );  
 iEvent.put(taus_pf_charge                               ,"tauspfcharge"                                   );  
 iEvent.put(taus_pf_sig_ncharge_cand                     ,"tauspfsignchargecand"                           );  
 iEvent.put(taus_pf_iso_ncharge_cand                     ,"tauspfisonchargecand"                           ); 
 iEvent.put(taus_pf_sig_nneutr_cand                      ,"tauspfsignneutrcand"                            ); 
 iEvent.put(taus_pf_iso_nneutr_cand                      ,"tauspfisonneutrcand"                            ); 
 iEvent.put(taus_pf_sig_ngamma_cand                      ,"tauspfsigngammacand"                            ); 
 iEvent.put(taus_pf_iso_ngamma_cand                      ,"tauspfisongammacand"                            ); 
 iEvent.put(taus_pf_lead_chargecand_p4                   ,"tauspfleadchargecandp4"                         ); 
 iEvent.put(taus_pf_lead_neutrcand_p4                    ,"tauspfleadneutrcandp4"                          ); 
 iEvent.put(taus_pf_lead_chargecand_Signed_Sipt          ,"tauspfleadchargecanSignedSipt"                  ); 
 iEvent.put(taus_pf_isolationchargecandPtSum             ,"tauspfisolationchargecandPtSum"                 );
 iEvent.put(taus_pf_isolationgammacandEtSum              ,"tauspfisolationgammacandEtSum"                  );
 iEvent.put(taus_pf_maximumHCALPFClusterEt               ,"tauspfmaximumHCALPFClusterEt"                   );
 iEvent.put(taus_pf_emf                                  ,"tauspfemf"                                      );
 iEvent.put(taus_pf_hcalTotOverPLead                     ,"tauspfhcalTotOverPLead"                         );
 iEvent.put(taus_pf_hcalMaxOverPLead                     ,"tauspfhcalMaxOverPLead"                         );
 iEvent.put(taus_pf_hcal3x3OverPLead                     ,"tauspfhcal3x3OverPLead"                         );
 
 
 iEvent.put(taus_pf_ecalStripSumEOverPLead               ,"tauspfecalStripSumEOverPLead"                   );
 iEvent.put(taus_pf_bremsRecoveryEOverPLead              ,"tauspfbremsRecoveryEOverPLead"                  );
 iEvent.put(taus_pf_electronPreID                        ,"tauspfelectronPreID"                            );
 
 iEvent.put(taus_pf_muonPreID                            ,"tauspfmuonPreID"                                );
 iEvent.put(taus_pf_hasMuonReference                     ,"tauspfhasMuonReference"                         );
 iEvent.put(taus_pf_caloComp                             ,"tauspfcaloComp"                                 );
 iEvent.put(taus_pf_segComp                              ,"tauspfsegComp"                                  );

 
 iEvent.put(taus_pf_leadtrk_p4                           ,"tauspfleadtrkp4"                                ); 
 iEvent.put(taus_pf_leadtrk_d0                           ,"tauspfleadtrkd0"                                );
 iEvent.put(taus_pf_leadtrk_z0                           ,"tauspfleadtrkz0"                                );
 iEvent.put(taus_pf_leadtrk_chi2                         ,"tauspfleadtrkchi2"                              );
 iEvent.put(taus_pf_leadtrk_ndof                         ,"tauspfleadtrkndof"                              );
 iEvent.put(taus_pf_leadtrk_validHits                    ,"tauspfleadtrkvalidHits"                         );
 iEvent.put(taus_pf_leadtrk_lostHits                     ,"tauspfleadtrklostHits"                          );
 
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFTauMaker);





  
