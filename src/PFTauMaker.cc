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
// $Id: PFTauMaker.cc,v 1.13 2010/04/25 13:56:53 kalavase Exp $
//
//


// system include files
#include <memory>
#include <math.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "CMS2/NtupleMaker/interface/PFTauMaker.h"
#include "CMS2/NtupleMaker/interface/CommonUtils.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;
using namespace CommonUtils;

//
// constructors and destructor
//

PFTauMaker::PFTauMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<vector<LorentzVector> >  (branchprefix+"p4"                            ).setBranchAlias(aliasprefix_+"_p4"                             );
  produces<vector<int> >            (branchprefix+"charge"                        ).setBranchAlias(aliasprefix_+"_charge"                         );

  produces<vector<vector <LorentzVector> > >  (branchprefix+"isochargecandp4"     ).setBranchAlias(aliasprefix_+"_isochargecand_p4"               );
  produces<vector<vector <LorentzVector> > >  (branchprefix+"isoneutrcandp4"      ).setBranchAlias(aliasprefix_+"_isoneutrcand_p4"                );
  produces<vector<vector <LorentzVector> > >  (branchprefix+"isogammacandp4"      ).setBranchAlias(aliasprefix_+"_isogammacand_p4"                );
  produces<vector<vector <LorentzVector> > >  (branchprefix+"sigchargecandp4"     ).setBranchAlias(aliasprefix_+"_sigchargecand_p4"               );
  produces<vector<vector <LorentzVector> > >  (branchprefix+"signeutrcandp4"      ).setBranchAlias(aliasprefix_+"_signeutrcand_p4"                );
  produces<vector<vector <LorentzVector> > >  (branchprefix+"siggammacandp4"      ).setBranchAlias(aliasprefix_+"_siggammacand_p4"                );

  

  produces<vector<LorentzVector> >  (branchprefix+"leadchargecandp4"              ).setBranchAlias(aliasprefix_+"_lead_chargecand_p4"             );
  produces<vector<LorentzVector> >  (branchprefix+"leadneutrcandp4"               ).setBranchAlias(aliasprefix_+"_lead_neutrcand_p4"              );
  produces<vector<float> >          (branchprefix+"leadchargecandSignedSipt"      ).setBranchAlias(aliasprefix_+"_lead_chargecand_Signed_Sipt"    );
  produces<vector<float> >          (branchprefix+"isolationchargecandPtSum"      ).setBranchAlias(aliasprefix_+"_isolationchargecandPtSum"       ); 
  produces<vector<float> >          (branchprefix+"isolationgammacandEtSum"       ).setBranchAlias(aliasprefix_+"_isolationgammacandEtSum"        ); 
  produces<vector<float> >          (branchprefix+"maximumHCALPFClusterEt"        ).setBranchAlias(aliasprefix_+"_maximumHCALPFClusterEt"         ); 
  //electron rejection
  produces<vector<float> >          (branchprefix+"emf"                           ).setBranchAlias(aliasprefix_+"_emf"                            ); 
  produces<vector<float> >          (branchprefix+"hcalTotOverPLead"              ).setBranchAlias(aliasprefix_+"_hcalTotOverPLead"               ); 
  produces<vector<float> >          (branchprefix+"hcalMaxOverPLead"              ).setBranchAlias(aliasprefix_+"_hcalMaxOverPLead"               ); 
  produces<vector<float> >          (branchprefix+"hcal3x3OverPLead"              ).setBranchAlias(aliasprefix_+"_hcal3x3OverPLead"               ); 
  produces<vector<float> >          (branchprefix+"ecalStripSumEOverPLead"        ).setBranchAlias(aliasprefix_+"_ecalStripSumEOverPLead"         ); 
  // produces<vector<float> >          (branchprefix+"bremsRecoveryEOverPLead"       ).setBranchAlias(aliasprefix_+"_bremsRecoveryEOverPLead"        ); 
  produces<vector<int> >            (branchprefix+"electronPreID"                 ).setBranchAlias(aliasprefix_+"_electronPreID"                  ); 
  produces<vector<float> >          (branchprefix+"electronPreIDOutput"           ).setBranchAlias(aliasprefix_+"_electronPreIDOutput"            ); 
  //muon rejection
  produces<vector<int> >            (branchprefix+"muonPreID"                     ).setBranchAlias(aliasprefix_+"_muonPreID"                      ); 
  produces<vector<int> >            (branchprefix+"hasMuonReference"              ).setBranchAlias(aliasprefix_+"_hasMuonReference"               ); 
  produces<vector<float> >          (branchprefix+"caloComp"                      ).setBranchAlias(aliasprefix_+"_caloComp"                       ); 
  produces<vector<float> >          (branchprefix+"segComp"                       ).setBranchAlias(aliasprefix_+"_segComp"                        ); 
  produces<vector<int> >            (branchprefix+"nmuonmatch"                    ).setBranchAlias(aliasprefix_+"_nmuonmatch"                     ); 

  //tau preID
  produces<vector<int> >            (branchprefix+"tightId"                       ).setBranchAlias(aliasprefix_+"_tightId"                        );

  produces<vector<int> >            (branchprefix+"leadtrkidx"                    ).setBranchAlias(aliasprefix_+"_leadtrk_idx"                    );
    
   
//get setup parameters
  pftausInputTag_              = iConfig.getParameter<InputTag>("pftausInputTag");
  minleadPFChargedHadrCandPt_  = iConfig.getParameter<double>("minleadPFChargedHadrCandPt");
  

}


PFTauMaker::~PFTauMaker() {}

void  PFTauMaker::beginJob() {
}

void PFTauMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void PFTauMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<LorentzVector> > taus_pf_p4                              (new vector<LorentzVector>            ) ;
  auto_ptr<vector<int> >           taus_pf_charge                          (new vector<int>                      ) ;

  auto_ptr<vector<LorentzVector> > taus_pf_lead_chargecand_p4              (new vector<LorentzVector>            ) ;
  auto_ptr<vector<LorentzVector> > taus_pf_lead_neutrcand_p4               (new vector<LorentzVector>            ) ;
  
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_isochargecand_p4       (new vector<vector<LorentzVector> >   ) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_isoneutrcand_p4        (new vector<vector<LorentzVector> >   ) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_isogammacand_p4        (new vector<vector<LorentzVector> >   ) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_sigchargecand_p4       (new vector<vector<LorentzVector> >   ) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_signeutrcand_p4        (new vector<vector<LorentzVector> >   ) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_siggammacand_p4        (new vector<vector<LorentzVector> >   ) ;
 
  auto_ptr<vector<float> >         taus_pf_lead_chargecand_Signed_Sipt     (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_isolationchargecandPtSum        (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_isolationgammacandEtSum         (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_maximumHCALPFClusterEt          (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_emf                             (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_hcalTotOverPLead                (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_hcalMaxOverPLead                (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_hcal3x3OverPLead                (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_ecalStripSumEOverPLead          (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_electronPreIDOutput             (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_caloComp                        (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_segComp                         (new vector<float>                    ) ;
  //auto_ptr<vector<float> >         taus_pf_bremsRecoveryEOverPLead         (new vector<float>) ;
 
  auto_ptr<vector<int> >           taus_pf_electronPreID                   (new vector<int>                      ) ;
  auto_ptr<vector<int> >           taus_pf_muonPreID                       (new vector<int>                      ) ;
  auto_ptr<vector<int> >           taus_pf_hasMuonReference                (new vector<int>                      ) ;
  auto_ptr<vector<int> >           taus_pf_nmuonmatch                      (new vector<int>                      ) ;
  auto_ptr<vector<int> >           taus_pf_tightId                         (new vector<int>                      ) ;
  auto_ptr<vector<int> >           taus_pf_leadtrk_idx                     (new vector<int>                      ) ;

 
  
  Handle<View<reco::PFTau> > taus_pf_h;
  iEvent.getByLabel(pftausInputTag_, taus_pf_h);
  
 
  size_t tausIndex = 0;
 for(View<reco::PFTau>::const_iterator tau_pf = taus_pf_h->begin();
      tau_pf != taus_pf_h->end(); tau_pf++, tausIndex++ ) {
   if(tau_pf->leadPFChargedHadrCand().isNull()) continue;
   if(tau_pf->leadPFChargedHadrCand()->pt()<minleadPFChargedHadrCandPt_) continue;

   taus_pf_p4                               ->push_back( LorentzVector( tau_pf->p4() ) );
   taus_pf_charge                           ->push_back( tau_pf->charge()              );

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

   if(pfIsoChargedCands.size() > 0){
     for(size_t iIsoCand = 0; iIsoCand < pfIsoChargedCands.size(); ++iIsoCand) 
       IsoChargedCands_p4.push_back( LorentzVector( pfIsoChargedCands[iIsoCand]->p4() ) );
   }

   if(pfIsoNeutrCands.size() >0) {
     for(size_t iIsoCand = 0; iIsoCand < pfIsoNeutrCands.size(); ++iIsoCand)
       IsoNeutrCands_p4.push_back( LorentzVector( pfIsoNeutrCands[iIsoCand]->p4() ) );
   }
   
   if(pfIsoGammaCands.size() >0) {
     for(size_t iIsoCand = 0; iIsoCand < pfIsoGammaCands.size(); ++iIsoCand) 
       IsoGammaCands_p4.push_back( LorentzVector( pfIsoGammaCands[iIsoCand]->p4() ) );
   }
   

   if(pfSigChargedCands.size() >0) {
     for(size_t iSigCand = 0; iSigCand < pfSigChargedCands.size(); ++iSigCand)
	 SigChargedCands_p4.push_back( LorentzVector( pfSigChargedCands[iSigCand]->p4() ) );
   }
   
  
   if(pfSigNeutrCands.size() >0){
     for(size_t iSigCand = 0; iSigCand < pfSigNeutrCands.size(); ++iSigCand) 
       SigNeutrCands_p4.push_back( LorentzVector( pfSigNeutrCands[iSigCand]->p4() ) );
   }

   
   if(pfSigGammaCands.size() >0){
     for(size_t iSigCand = 0; iSigCand < pfSigGammaCands.size(); ++iSigCand)
       SigGammaCands_p4.push_back( LorentzVector( pfSigGammaCands[iSigCand]->p4() ) );	 
   }
   
   
   
   taus_pf_isochargecand_p4 -> push_back( IsoChargedCands_p4 );
   taus_pf_isoneutrcand_p4  -> push_back( IsoNeutrCands_p4   );
   taus_pf_isogammacand_p4  -> push_back( IsoGammaCands_p4   );
   taus_pf_sigchargecand_p4 -> push_back( SigChargedCands_p4 );
   taus_pf_signeutrcand_p4  -> push_back( SigNeutrCands_p4   );
   taus_pf_siggammacand_p4  -> push_back( SigGammaCands_p4   );
   
   taus_pf_lead_chargecand_p4               ->push_back( LorentzVector( tau_pf->leadPFChargedHadrCand().get()->p4() ) );
   
   taus_pf_lead_neutrcand_p4                ->push_back( tau_pf->leadPFNeutralCand().isNonnull()? LorentzVector( tau_pf->leadPFNeutralCand().get()->p4() ) :  LorentzVector(0, 0, 0, 0)  );
   taus_pf_lead_chargecand_Signed_Sipt      ->push_back( !isfinite(tau_pf->leadPFChargedHadrCandsignedSipt()) ? 
							 -9999 : tau_pf->leadPFChargedHadrCandsignedSipt()  ); 
   taus_pf_isolationchargecandPtSum         ->push_back( !isfinite(tau_pf->isolationPFChargedHadrCandsPtSum()) ? 
							 -9999 : tau_pf->isolationPFChargedHadrCandsPtSum() ); 
   taus_pf_isolationgammacandEtSum          ->push_back( !isfinite(tau_pf->isolationPFGammaCandsEtSum())       ? 
							 -9999. : tau_pf->isolationPFGammaCandsEtSum()      ); 
   taus_pf_maximumHCALPFClusterEt           ->push_back( !isfinite(tau_pf->maximumHCALPFClusterEt())           ? 
							 -9999. : tau_pf->maximumHCALPFClusterEt()          ); 
   taus_pf_emf                              ->push_back( !isfinite(tau_pf->emFraction())                       ? 
							 -9999. : tau_pf->emFraction()                      ); 
   taus_pf_hcalTotOverPLead                 ->push_back( !isfinite(tau_pf->hcalTotOverPLead())                  ?
							 -9999. : tau_pf->hcalTotOverPLead()                ); 
   taus_pf_hcalMaxOverPLead                 ->push_back( !isfinite(tau_pf->hcalMaxOverPLead())                 ?
							 -9999. : tau_pf->hcalMaxOverPLead()                ); 
   taus_pf_hcal3x3OverPLead                 ->push_back( !isfinite(tau_pf->hcal3x3OverPLead())                 ?
							 -9999. : tau_pf->hcal3x3OverPLead()                );
   taus_pf_ecalStripSumEOverPLead           ->push_back( !isfinite(tau_pf->ecalStripSumEOverPLead())           ? 
							 -9999. : tau_pf->ecalStripSumEOverPLead()          ); 
   //taus_pf_bremsRecoveryEOverPLead          ->push_back( tau_pf->bremsRecoveryEOverPLead()          ); 
 
  
   if(tau_pf->electronPreIDDecision()) taus_pf_electronPreID      ->push_back(1);
   else                                taus_pf_electronPreID      ->push_back(0);
   if(!tau_pf->muonDecision())         taus_pf_muonPreID          ->push_back(1);
   else                                taus_pf_muonPreID          ->push_back(0);
  
   
   taus_pf_electronPreIDOutput      ->push_back(tau_pf->electronPreIDOutput());
   if(tau_pf->hasMuonReference()){
                                      taus_pf_hasMuonReference   ->push_back(1);
                                      taus_pf_caloComp           ->push_back(!isfinite(tau_pf->caloComp()) ? 
									     -9999. : tau_pf->caloComp()  );
                                      taus_pf_segComp            ->push_back(!isfinite(tau_pf->segComp()) ? 
									     -9999. : tau_pf->segComp() );
				      MuonRef muonref = tau_pf->leadPFChargedHadrCand()->muonRef();
				      taus_pf_nmuonmatch      ->push_back(!isfinite(muonref ->numberOfMatches()) ?
									  -9999. : muonref->numberOfMatches());
   }
   else {
                                      taus_pf_hasMuonReference   ->push_back(0);
                                      taus_pf_caloComp           ->push_back(-9999.);
                                      taus_pf_segComp            ->push_back(-9999.);
				      taus_pf_nmuonmatch         ->push_back(0);
				     
   }
   
   const edm::RefToBase<reco::PFTau> pftauRef =  taus_pf_h->refAt(tausIndex);
   bool tight_id = identify(pftauRef) ;
   if(tight_id) taus_pf_tightId                         ->push_back(1);
   else         taus_pf_tightId                         ->push_back(0);
   
   const TrackRef leadTrack = tau_pf->leadTrack()  ;
   if(leadTrack.isNonnull()){
     taus_pf_leadtrk_idx                    ->push_back( static_cast<int>(leadTrack.key())                      );
   }
   else {
     taus_pf_leadtrk_idx                    ->push_back( -9999                                               );
      
   }
  
   }

 std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");


 iEvent.put(taus_pf_isochargecand_p4                    ,branchprefix+"isochargecandp4"                           );  
 iEvent.put(taus_pf_isoneutrcand_p4                     ,branchprefix+"isoneutrcandp4"                            );  
 iEvent.put(taus_pf_isogammacand_p4                     ,branchprefix+"isogammacandp4"                            );  
 iEvent.put(taus_pf_sigchargecand_p4                    ,branchprefix+"sigchargecandp4"                           );  
 iEvent.put(taus_pf_signeutrcand_p4                     ,branchprefix+"signeutrcandp4"                            );  
 iEvent.put(taus_pf_siggammacand_p4                     ,branchprefix+"siggammacandp4"                            );  
 
 iEvent.put(taus_pf_p4                                   ,branchprefix+"p4"                                       );  
 iEvent.put(taus_pf_charge                               ,branchprefix+"charge"                                   );  
 iEvent.put(taus_pf_lead_chargecand_p4                   ,branchprefix+"leadchargecandp4"                         ); 
 iEvent.put(taus_pf_lead_neutrcand_p4                    ,branchprefix+"leadneutrcandp4"                          ); 
 iEvent.put(taus_pf_lead_chargecand_Signed_Sipt          ,branchprefix+"leadchargecandSignedSipt"                 ); 
 iEvent.put(taus_pf_isolationchargecandPtSum             ,branchprefix+"isolationchargecandPtSum"                 );
 iEvent.put(taus_pf_isolationgammacandEtSum              ,branchprefix+"isolationgammacandEtSum"                  );
 iEvent.put(taus_pf_maximumHCALPFClusterEt               ,branchprefix+"maximumHCALPFClusterEt"                   );
 iEvent.put(taus_pf_emf                                  ,branchprefix+"emf"                                      );
 iEvent.put(taus_pf_hcalTotOverPLead                     ,branchprefix+"hcalTotOverPLead"                         );
 iEvent.put(taus_pf_hcalMaxOverPLead                     ,branchprefix+"hcalMaxOverPLead"                         );
 iEvent.put(taus_pf_hcal3x3OverPLead                     ,branchprefix+"hcal3x3OverPLead"                         );
 
 
 iEvent.put(taus_pf_ecalStripSumEOverPLead               ,branchprefix+"ecalStripSumEOverPLead"                   );
 //iEvent.put(taus_pf_bremsRecoveryEOverPLead              ,branchprefix+"bremsRecoveryEOverPLead"                  );
 iEvent.put(taus_pf_electronPreID                        ,branchprefix+"electronPreID"                            );
 iEvent.put(taus_pf_electronPreIDOutput                  ,branchprefix+"electronPreIDOutput"                      );
 
 iEvent.put(taus_pf_muonPreID                            ,branchprefix+"muonPreID"                                );
 iEvent.put(taus_pf_hasMuonReference                     ,branchprefix+"hasMuonReference"                         );
 iEvent.put(taus_pf_caloComp                             ,branchprefix+"caloComp"                                 );
 iEvent.put(taus_pf_segComp                              ,branchprefix+"segComp"                                  );
 iEvent.put(taus_pf_nmuonmatch                           ,branchprefix+"nmuonmatch"                               );

 iEvent.put(taus_pf_tightId                              ,branchprefix+"tightId"                                  );
 iEvent.put(taus_pf_leadtrk_idx                          ,branchprefix+"leadtrkidx"                               ); 


  
 
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





  
