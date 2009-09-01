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
// $Id: PFTauMaker.cc,v 1.3 2009/09/01 00:44:20 yanjuntu Exp $
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
#include "DataFormats/MuonReco/interface/Muon.h"

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

//   produces<vector<int> >            ("tauspfsignchargecand"                ).setBranchAlias("taus_pf_sig_ncharge_cand"               );
//   produces<vector<int> >            ("tauspfisonchargecand"                ).setBranchAlias("taus_pf_iso_ncharge_cand"               );
//   produces<vector<int> >            ("tauspfsignneutrcand"                 ).setBranchAlias("taus_pf_sig_nneutr_cand"                );
//   produces<vector<int> >            ("tauspfisonneutrcand"                 ).setBranchAlias("taus_pf_iso_nneutr_cand"                );
//   produces<vector<int> >            ("tauspfsigngammacand"                 ).setBranchAlias("taus_pf_sig_ngamma_cand"                );
//   produces<vector<int> >            ("tauspfisongammacand"                 ).setBranchAlias("taus_pf_iso_ngamma_cand"                );

  produces<vector<vector <LorentzVector> > >  ("tauspfisochargecandp4"     ).setBranchAlias("taus_pf_isochargecand_p4"               );
  produces<vector<vector <LorentzVector> > >  ("tauspfisoneutrcandp4"      ).setBranchAlias("taus_pf_isoneutr_p4"                    );
  produces<vector<vector <LorentzVector> > >  ("tauspfisogammacandp4"      ).setBranchAlias("taus_pf_isogammacand_p4"                );
  produces<vector<vector <LorentzVector> > >  ("tauspfsigchargecandp4"     ).setBranchAlias("taus_pf_sigchargecand_p4"               );
  produces<vector<vector <LorentzVector> > >  ("tauspfsigneutrcandp4"      ).setBranchAlias("taus_pf_signeutr_p4"                    );
  produces<vector<vector <LorentzVector> > >  ("tauspfsiggammacandp4"      ).setBranchAlias("taus_pf_siggammacand_p4"                );

  

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
  produces<vector<float> >          ("tauspfelectronPreIDOutput"           ).setBranchAlias("taus_pf_electronPreIDOutput"            ); 
  //muon rejection
  produces<vector<int> >            ("tauspfmuonPreID"                     ).setBranchAlias("taus_pf_muonPreID"                      ); 
  produces<vector<int> >            ("tauspfhasMuonReference"              ).setBranchAlias("taus_pf_hasMuonReference"               ); 
  produces<vector<float> >          ("tauspfcaloComp"                      ).setBranchAlias("taus_pf_caloComp"                       ); 
  produces<vector<float> >          ("tauspfsegComp"                       ).setBranchAlias("taus_pf_segComp"                        ); 
  produces<vector<int> >            ("tauspfnmuonmatch"                    ).setBranchAlias("taus_pf_nmuonmatch"                     ); 

  //tau preID
  produces<vector<int> >            ("tauspftightId"                       ).setBranchAlias("taus_pf_tightId"                        );

  
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

//   auto_ptr<vector<int> >           taus_pf_sig_ncharge_cand                (new vector<int>) ;
//   auto_ptr<vector<int> >           taus_pf_iso_ncharge_cand                (new vector<int>) ;
//   auto_ptr<vector<int> >           taus_pf_sig_nneutr_cand                 (new vector<int>) ;
//   auto_ptr<vector<int> >           taus_pf_iso_nneutr_cand                 (new vector<int>) ;
//   auto_ptr<vector<int> >           taus_pf_sig_ngamma_cand                 (new vector<int>) ;
//   auto_ptr<vector<int> >           taus_pf_iso_ngamma_cand                 (new vector<int>) ;
  auto_ptr<vector<LorentzVector> > taus_pf_lead_chargecand_p4              (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > taus_pf_lead_neutrcand_p4               (new vector<LorentzVector>) ;
  
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_isochargecand_p4       (new vector<vector<LorentzVector> >) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_isoneutrcand_p4        (new vector<vector<LorentzVector> >) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_isogammacand_p4        (new vector<vector<LorentzVector> >) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_sigchargecand_p4       (new vector<vector<LorentzVector> >) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_signeutrcand_p4        (new vector<vector<LorentzVector> >) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_siggammacand_p4        (new vector<vector<LorentzVector> >) ;
 
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
  auto_ptr<vector<float> >         taus_pf_electronPreIDOutput             (new vector<float>) ;
  auto_ptr<vector<int> >           taus_pf_muonPreID                       (new vector<int>) ;
  auto_ptr<vector<int> >           taus_pf_hasMuonReference                (new vector<int>) ;
  auto_ptr<vector<float> >         taus_pf_caloComp                        (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_segComp                         (new vector<float>) ;
  auto_ptr<vector<int> >           taus_pf_nmuonmatch                      (new vector<int>) ;

  auto_ptr<vector<int> >           taus_pf_tightId                         (new vector<int>) ;
 
  auto_ptr<vector<LorentzVector> > taus_pf_leadtrk_p4                      (new vector<LorentzVector>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_d0                      (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_z0                      (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_chi2                    (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_ndof                    (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_validHits               (new vector<float>) ;
  auto_ptr<vector<float> >         taus_pf_leadtrk_lostHits                (new vector<float>) ;
 

 
  
  Handle<View<reco::PFTau> > taus_pf_h;
  iEvent.getByLabel(pftausInputTag, taus_pf_h);
  
 
  size_t tausIndex = 0;
 for(View<reco::PFTau>::const_iterator tau_pf = taus_pf_h->begin();
      tau_pf != taus_pf_h->end(); tau_pf++, tausIndex++ ) {
   if(tau_pf->leadPFChargedHadrCand().isNull()) continue;
   if(tau_pf->leadPFChargedHadrCand()->pt()<5.0) continue;

   //printf("%s  \n", "PFTau  ");

   taus_pf_p4                               ->push_back( tau_pf->p4()                                       );
   taus_pf_charge                           ->push_back( tau_pf->charge()                                   );
//    taus_pf_sig_ncharge_cand                 ->push_back( tau_pf->signalPFChargedHadrCands().size()          ); 
//    taus_pf_iso_ncharge_cand                 ->push_back( tau_pf->isolationPFChargedHadrCands().size()       );                                                                            
//    taus_pf_sig_nneutr_cand                  ->push_back( tau_pf->signalPFNeutrHadrCands().size()            ); 
//    taus_pf_iso_nneutr_cand                  ->push_back( tau_pf->isolationPFNeutrHadrCands().size()         ); 
//    taus_pf_sig_ngamma_cand                  ->push_back( tau_pf->signalPFGammaCands().size()                ); 
//    taus_pf_iso_ngamma_cand                  ->push_back( tau_pf->isolationPFGammaCands().size()             ); 

   vector<LorentzVector>  IsoChargedCands_p4;
   vector<LorentzVector>  IsoNeutrCands_p4;
   vector<LorentzVector>  IsoGammaCands_p4;
   vector<LorentzVector>  SigChargedCands_p4;
   vector<LorentzVector>  SigNeutrCands_p4;
   vector<LorentzVector>  SigGammaCands_p4;
 
   const PFCandidateRefVector& pfIsoChargedCands =  tau_pf->isolationPFChargedHadrCands();
   const PFCandidateRefVector& pfIsoNeutrCands   =  tau_pf->isolationPFNeutrHadrCands();
   const PFCandidateRefVector& pfIsoGammaCands   =  tau_pf->isolationPFGammaCands();
   const PFCandidateRefVector& pfSigChargedCands =  tau_pf->signalPFChargedHadrCands();
   const PFCandidateRefVector& pfSigNeutrCands   =  tau_pf->signalPFNeutrHadrCands();
   const PFCandidateRefVector& pfSigGammaCands   =  tau_pf->signalPFGammaCands();

   if(pfIsoChargedCands.size() >0){
     for(size_t iIsoCand = 0; iIsoCand < pfIsoChargedCands.size(); ++iIsoCand)
       {
	 
	 IsoChargedCands_p4.push_back( pfIsoChargedCands[iIsoCand]->p4()) ;
	 
       }
   }
   else IsoChargedCands_p4.push_back( LorentzVector(0, 0, 0, 0) );

  
   if(pfIsoNeutrCands.size() >0){
     for(size_t iIsoCand = 0; iIsoCand < pfIsoNeutrCands.size(); ++iIsoCand)
       {
	 
	 IsoNeutrCands_p4.push_back( pfIsoNeutrCands[iIsoCand]->p4()) ;
	 
       }
   }
   else IsoNeutrCands_p4.push_back( LorentzVector(0, 0, 0, 0) );
   
   if(pfIsoGammaCands.size() >0){
     for(size_t iIsoCand = 0; iIsoCand < pfIsoGammaCands.size(); ++iIsoCand)
       {
	 
	 IsoGammaCands_p4.push_back( pfIsoGammaCands[iIsoCand]->p4()) ;
	 
       }
   }
   else IsoGammaCands_p4.push_back( LorentzVector(0, 0, 0, 0) );
   

  if(pfSigChargedCands.size() >0){
     for(size_t iSigCand = 0; iSigCand < pfSigChargedCands.size(); ++iSigCand)
       {
	 
	 SigChargedCands_p4.push_back( pfSigChargedCands[iSigCand]->p4()) ;
	 
       }
   }
   else SigChargedCands_p4.push_back( LorentzVector(0, 0, 0, 0) );

  
   if(pfSigNeutrCands.size() >0){
     for(size_t iSigCand = 0; iSigCand < pfSigNeutrCands.size(); ++iSigCand)
       {
	 
	 SigNeutrCands_p4.push_back( pfSigNeutrCands[iSigCand]->p4()) ;
	 
       }
   }
   else SigNeutrCands_p4.push_back( LorentzVector(0, 0, 0, 0) );
   
   if(pfSigGammaCands.size() >0){
     for(size_t iSigCand = 0; iSigCand < pfSigGammaCands.size(); ++iSigCand)
       {
	 
	 SigGammaCands_p4.push_back( pfSigGammaCands[iSigCand]->p4()) ;
	 
       }
   }
   else SigGammaCands_p4.push_back( LorentzVector(0, 0, 0, 0) );
   
   

   taus_pf_isochargecand_p4 -> push_back( IsoChargedCands_p4 );
   taus_pf_isoneutrcand_p4 -> push_back( IsoNeutrCands_p4 );
   taus_pf_isogammacand_p4 -> push_back( IsoGammaCands_p4 );
   taus_pf_sigchargecand_p4 -> push_back( SigChargedCands_p4 );
   taus_pf_signeutrcand_p4 -> push_back( SigNeutrCands_p4 );
   taus_pf_siggammacand_p4 -> push_back( SigGammaCands_p4 );
   
   taus_pf_lead_chargecand_p4               ->push_back( tau_pf->leadPFChargedHadrCand().get()->p4()        );
   
   taus_pf_lead_neutrcand_p4                ->push_back( tau_pf->leadPFNeutralCand().isNonnull()? tau_pf->leadPFNeutralCand().get()->p4() :  LorentzVector(0, 0, 0, 0)  );
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
  
   taus_pf_electronPreIDOutput      ->push_back(tau_pf->electronPreIDOutput());
   if(tau_pf->hasMuonReference()){
                                      taus_pf_hasMuonReference   ->push_back(1);
                                      taus_pf_caloComp           ->push_back(tau_pf->caloComp());
                                      taus_pf_segComp            ->push_back(tau_pf->segComp());
				      MuonRef muonref = tau_pf->leadPFChargedHadrCand()->muonRef();
				      taus_pf_nmuonmatch         ->push_back(muonref ->numberOfMatches());
   }
   else {
                                      taus_pf_hasMuonReference   ->push_back(0);
                                      taus_pf_caloComp           ->push_back(-999.);
                                      taus_pf_segComp            ->push_back(-999.);
				      taus_pf_nmuonmatch         ->push_back(0);
				     
   }
   
   const edm::RefToBase<reco::PFTau> pftauRef =  taus_pf_h->refAt(tausIndex);
   bool tight_id = identify(pftauRef) ;
   if(tight_id) taus_pf_tightId                         ->push_back(1);
   else         taus_pf_tightId                         ->push_back(0);
   
   const TrackRef leadTrack = tau_pf->leadTrack()  ;
   if(leadTrack.isNonnull()){
     taus_pf_leadtrk_p4                    ->push_back( LorentzVector( leadTrack.get()->px(), leadTrack.get()->py(),
								       leadTrack.get()->pz(), leadTrack.get()->p())                );
     taus_pf_leadtrk_d0                    ->push_back( leadTrack->d0()                                       );
     taus_pf_leadtrk_z0                    ->push_back( leadTrack->dz()                                       );
     taus_pf_leadtrk_chi2                  ->push_back( leadTrack->chi2()                                     );
     taus_pf_leadtrk_ndof                  ->push_back( leadTrack->ndof()                                     );
     taus_pf_leadtrk_validHits             ->push_back( leadTrack->numberOfValidHits()                        );
     taus_pf_leadtrk_lostHits              ->push_back( leadTrack->numberOfLostHits()                         );
   }
   else {
     taus_pf_leadtrk_p4                    ->push_back( LorentzVector( 0, 0,0,0 )                             );
     taus_pf_leadtrk_d0                    ->push_back( -999.                                                 );
     taus_pf_leadtrk_z0                    ->push_back( -999.                                                 );
     taus_pf_leadtrk_chi2                  ->push_back( -999.                                                 );
     taus_pf_leadtrk_ndof                  ->push_back( -999.                                                 );
     taus_pf_leadtrk_validHits             ->push_back( -999.                                                 );
     taus_pf_leadtrk_lostHits              ->push_back( -999.                                                 );
     
   }
  
 }

 iEvent.put(taus_pf_isochargecand_p4                    ,"tauspfisochargecandp4"                           );  
 iEvent.put(taus_pf_isoneutrcand_p4                     ,"tauspfisoneutrcandp4"                            );  
 iEvent.put(taus_pf_isogammacand_p4                     ,"tauspfisogammacandp4"                            );  
 iEvent.put(taus_pf_sigchargecand_p4                    ,"tauspfsigchargecandp4"                           );  
 iEvent.put(taus_pf_signeutrcand_p4                     ,"tauspfsigneutrcandp4"                            );  
 iEvent.put(taus_pf_siggammacand_p4                     ,"tauspfsiggammacandp4"                            );  
 
 iEvent.put(taus_pf_p4                                   ,"tauspfp4"                                       );  
 iEvent.put(taus_pf_charge                               ,"tauspfcharge"                                   );  
//  iEvent.put(taus_pf_sig_ncharge_cand                     ,"tauspfsignchargecand"                           );  
//  iEvent.put(taus_pf_iso_ncharge_cand                     ,"tauspfisonchargecand"                           ); 
//  iEvent.put(taus_pf_sig_nneutr_cand                      ,"tauspfsignneutrcand"                            ); 
//  iEvent.put(taus_pf_iso_nneutr_cand                      ,"tauspfisonneutrcand"                            ); 
//  iEvent.put(taus_pf_sig_ngamma_cand                      ,"tauspfsigngammacand"                            ); 
//  iEvent.put(taus_pf_iso_ngamma_cand                      ,"tauspfisongammacand"                            ); 
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
 iEvent.put(taus_pf_electronPreIDOutput                  ,"tauspfelectronPreIDOutput"                            );
 
 iEvent.put(taus_pf_muonPreID                            ,"tauspfmuonPreID"                                );
 iEvent.put(taus_pf_hasMuonReference                     ,"tauspfhasMuonReference"                         );
 iEvent.put(taus_pf_caloComp                             ,"tauspfcaloComp"                                 );
 iEvent.put(taus_pf_segComp                              ,"tauspfsegComp"                                  );
 iEvent.put(taus_pf_nmuonmatch                            ,"tauspfnmuonmatch"                                );

 iEvent.put(taus_pf_tightId                              ,"tauspftightId"                                  );
 

 

 iEvent.put(taus_pf_leadtrk_p4                           ,"tauspfleadtrkp4"                                ); 
 iEvent.put(taus_pf_leadtrk_d0                           ,"tauspfleadtrkd0"                                );
 iEvent.put(taus_pf_leadtrk_z0                           ,"tauspfleadtrkz0"                                );
 iEvent.put(taus_pf_leadtrk_chi2                         ,"tauspfleadtrkchi2"                              );
 iEvent.put(taus_pf_leadtrk_ndof                         ,"tauspfleadtrkndof"                              );
 iEvent.put(taus_pf_leadtrk_validHits                    ,"tauspfleadtrkvalidHits"                         );
 iEvent.put(taus_pf_leadtrk_lostHits                     ,"tauspfleadtrklostHits"                          );
 
 
}

//-------------------------------------------------------------------------------------------------
//Tau ID
//-------------------------------------------------------------------------------------------------
bool PFTauMaker::identify(const edm::RefToBase<reco::PFTau> &tau_pf) {
  /*
  float emFraction_maxValue_ = 0.9;
  float hcalTotOverPLead_minValue_ = 0.1;
  float hcal3x3OverPLead_minValue_ = 0.1;
  float hcalMaxOverPLead_minValue_ = 0.1;
  float  EOverPLead_maxValue_ = 1.8;
  float  EOverPLead_minValue_ = 0.8; 
  float  bremsRecoveryEOverPLead_minValue_ = 0.8;
  float  bremsRecoveryEOverPLead_maxValue_ = 1.8; 
  float  elecPreID0_EOverPLead_maxValue_ = 0.95;
  float  elecPreID0_HOverPLead_minValue_ = 0.05;
  float  elecPreID1_EOverPLead_maxValue_ = 0.8;
  float  elecPreID1_HOverPLead_minValue_ = 0.15;
  */
  float  pfelectronMVA_maxValue_ = -0.1;
  
  //isolation
  
  const PFCandidateRefVector& pfIsoChargedCands =  tau_pf->isolationPFChargedHadrCands();
  unsigned int tracksAboveThreshold = 0;
  for(size_t iIsoCand = 0; iIsoCand < pfIsoChargedCands.size(); ++iIsoCand)
    {
      if(pfIsoChargedCands[iIsoCand]->pt() > 1.0) {
	if(++tracksAboveThreshold > 0) {
	  return false;
	
	}
      }
    }
  /*
  const TrackRefVector& isolationTracks = tau_pf->isolationTracks();
  unsigned int tracksAboveThreshold = 0;
  for(size_t iTrack = 0; iTrack < isolationTracks.size(); ++iTrack)
    {
      if(isolationTracks[iTrack]->pt() > 1.0) {
	if(++tracksAboveThreshold > 0)
	  {
	 
	    return false;
	  }
      }
    }
  */
  const PFCandidateRefVector& pfIsoGammaCands = tau_pf->isolationPFGammaCands();
  unsigned int gammasAboveThreshold = 0;
  for(size_t iIsoGamma = 0; iIsoGamma < pfIsoGammaCands.size(); ++iIsoGamma)
    {
      if(pfIsoGammaCands[iIsoGamma]->pt() > 1.5) {
	if(++gammasAboveThreshold > 0) {
	  return false;
	}
      }
    }	 

  
  
  //anti-electron
  bool decision = false;
  bool emfPass = true, htotPass = true, hmaxPass = true; 
  bool h3x3Pass = true, estripPass = true, erecovPass = true;
  bool epreidPass = true, epreid2DPass = true;
  bool mvaPass = true;

  //optional ways to reject electrons from taus. For the time being, we use the default



  /*
  //cout<<"emfraction  = "<<tau_pf->emFraction()<<endl;
  if ( tau_pf->emFraction() > emFraction_maxValue_) {
     emfPass = false;
   }

  if (tau_pf->hcalTotOverPLead() < hcalTotOverPLead_minValue_) {
    htotPass = false;
  }

  if (tau_pf->hcalMaxOverPLead() < hcalMaxOverPLead_minValue_) {
    hmaxPass = false;
  }

  if (tau_pf->hcal3x3OverPLead() < hcal3x3OverPLead_minValue_) {
    h3x3Pass = false;
  }

  if (tau_pf->ecalStripSumEOverPLead() > EOverPLead_minValue_ &&
      tau_pf->ecalStripSumEOverPLead() < EOverPLead_maxValue_) {
    estripPass = false;
  } else {
    estripPass = true;
  }

    if (tau_pf->bremsRecoveryEOverPLead() > bremsRecoveryEOverPLead_minValue_ &&
	tau_pf->bremsRecoveryEOverPLead() < bremsRecoveryEOverPLead_maxValue_) {
      erecovPass = false;
    } 
    else {
      erecovPass = true;
    } 

    if (tau_pf->electronPreIDDecision()) {
      epreidPass = false;
    }  
    else {
      epreidPass = true;
    }

    if (
	(tau_pf->electronPreIDDecision()  && (tau_pf->ecalStripSumEOverPLead() < elecPreID1_EOverPLead_maxValue_ ||tau_pf->hcal3x3OverPLead() > elecPreID1_HOverPLead_minValue_)) ||
	(!tau_pf->electronPreIDDecision() && (tau_pf->ecalStripSumEOverPLead() < elecPreID0_EOverPLead_maxValue_ ||tau_pf->hcal3x3OverPLead() > elecPreID0_HOverPLead_minValue_))
	)
      {
	epreid2DPass = true;
      }  else {
	epreid2DPass = false;
      }
  */
    if (tau_pf->electronPreIDOutput()>pfelectronMVA_maxValue_) {
       mvaPass = false;
     
     }  
  
    decision = emfPass && htotPass && hmaxPass && 
      h3x3Pass && estripPass && erecovPass && epreidPass && epreid2DPass && mvaPass;

    
    if(!decision) return false;
  
    //anti-muon

    bool decision2 = true;
 
    if(tau_pf->hasMuonReference() ){
    
      MuonRef muonref = tau_pf->leadPFChargedHadrCand()->muonRef();
      if ( muonref ->numberOfMatches() > 0 ) {
	decision2 = false;
      }
    }
    
    if(!decision2) return false;

    return true;
}
//define this as a plug-in
DEFINE_FWK_MODULE(PFTauMaker);





  
