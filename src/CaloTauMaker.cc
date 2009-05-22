//-*- C++ -*-
//
// Package:    CaloTauMaker
// Class:      CaloTauMaker
// 
/**\class CaloTauMaker CaloTauMaker.cc CMS2/NtupleMaker/src/CaloTauMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// $Id: CaloTauMaker.cc,v 1.1 2009/05/22 18:34:31 yanjuntu Exp $
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

#include "CMS2/NtupleMaker/interface/CaloTauMaker.h"
#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/CaloTauFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"


typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

CaloTauMaker::CaloTauMaker(const edm::ParameterSet& iConfig) {



  produces<vector<LorentzVector> >  ("tauscalop4"                            ).setBranchAlias("taus_calo_p4"                      );
  produces<vector<int> >            ("tauscalosigntrks"                      ).setBranchAlias("taus_calo_sig_ntrks"               );
  produces<vector<int> >            ("tauscaloisontrks"                      ).setBranchAlias("taus_calo_iso_ntrks"               );
 
  produces<vector<LorentzVector> >  ("tauscaloleadtrkp4"                     ).setBranchAlias("taus_calo_leadtrk_p4"              );
  produces<vector<int> >            ("tauscalocharge"                        ).setBranchAlias("taus_calo_charge"                  );
  produces<vector<float> >          ("tauscaloleadtrkd0"                     ).setBranchAlias("taus_calo_leadtrk_d0"              );  
  produces<vector<float> >          ("tauscaloleadtrkz0"                     ).setBranchAlias("taus_calo_leadtrk_z0"              ); 
  produces<vector<float> >          ("tauscaloleadtrkchi2"                   ).setBranchAlias("taus_calo_leadtrk_chi2"            ); 
  produces<vector<float> >          ("tauscaloleadtrkndof"                   ).setBranchAlias("taus_calo_leadtrk_ndof"            ); 
  produces<vector<float> >          ("tauscaloleadtrkvalidHits"              ).setBranchAlias("taus_calo_leadtrk_validHits"       ); 
  produces<vector<float> >          ("tauscaloleadtrklostHits"               ).setBranchAlias("taus_calo_leadtrk_lostHits"        ); 
  produces<vector<float> >          ("tauscaloleadtrkSignedSipt"             ).setBranchAlias("taus_calo_leadtrk_Signed_Sipt"     );  
  produces<vector<float> >          ("tauscaloleadtrkHCAL3x3hitsEtSum"       ).setBranchAlias("taus_calo_leadtrk_HCAL3x3hitsEtSum"); 
  produces<vector<float> >          ("tauscaloleadtrkHCAL3x3hottesthitDEta"  ).setBranchAlias("taus_calo_leadtrk_HCAL3x3hottesthitDEta"); 
  
  produces<vector<float> >          ("tauscalosignaltrksInvariantMass"        ).setBranchAlias("taus_calo_signaltrksInvariantMass" ); 
  produces<vector<float> >          ("tauscaloisolationtrksPtSum"             ).setBranchAlias("taus_calo_isolationtrksPtSum"      ); 
  produces<vector<float> >          ("tauscaloisolationECALhitsEtSum"         ).setBranchAlias("taus_calo_isolationECALhitsEtSum"  ); 
  produces<vector<float> >          ("tauscalomaximumHCALhitEt"               ).setBranchAlias("taus_calo_maximumHCALhitEt"        ); 
  
  

  //leadTrackECALSurfContactPoint, leadTrackavoidsECALcrack
  
  
    
   
//get setup parameters

  calotausInputTag    = iConfig.getParameter<InputTag>("calotausInputTag");
}


CaloTauMaker::~CaloTauMaker() {}

void  CaloTauMaker::beginJob(const edm::EventSetup&) {
}

void CaloTauMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void CaloTauMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
 
  
  auto_ptr<vector<LorentzVector> > taus_calo_p4                              (new vector<LorentzVector>) ;
  auto_ptr<vector<int> >           taus_calo_sig_ntrks                       (new vector<int>) ;
  auto_ptr<vector<int> >           taus_calo_iso_ntrks                       (new vector<int>) ;
  auto_ptr<vector<int> >           taus_calo_charge                          (new vector<int>) ;
  auto_ptr<vector<LorentzVector> > taus_calo_leadtrk_p4                      (new vector<LorentzVector>) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_d0                      (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_z0                      (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_chi2                    (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_ndof                    (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_validHits               (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_lostHits                (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_Signed_Sipt             (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_HCAL3x3hitsEtSum        (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_HCAL3x3hottesthitDEta   (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_signaltrksInvariantMass         (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_isolationtrksPtSum              (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_isolationECALhitsEtSum          (new vector<float>) ;
  auto_ptr<vector<float> >         taus_calo_maximumHCALhitEt                (new vector<float>) ;
  
  
  
  Handle<View<reco::CaloTau> > taus_calo_h;
  iEvent.getByLabel(calotausInputTag, taus_calo_h);

  
 

 for(View<reco::CaloTau>::const_iterator tau_calo = taus_calo_h->begin();
      tau_calo != taus_calo_h->end(); tau_calo++) {
   if(tau_calo->leadTrack().isNull())  continue;
   if(tau_calo->leadTrack()->pt()<5.0) continue;
   //printf("%s  \n", "CaloTau  ");
   
   const TrackRef leadTrack = tau_calo->leadTrack()  ;
   
   taus_calo_p4                            ->push_back( tau_calo->p4()                                        );
   taus_calo_sig_ntrks                     ->push_back( tau_calo->signalTracks().size()                       );
   taus_calo_iso_ntrks                     ->push_back( tau_calo->isolationTracks().size()                    );
   taus_calo_leadtrk_p4                    ->push_back( LorentzVector( leadTrack.get()->px(), leadTrack.get()->py(),
						leadTrack.get()->pz(), leadTrack.get()->p())                  );
   taus_calo_charge                        ->push_back( tau_calo->charge()                                    );
   taus_calo_leadtrk_d0                    ->push_back( leadTrack->d0()                                       );
   taus_calo_leadtrk_z0                    ->push_back( leadTrack->dz()                                       );
   taus_calo_leadtrk_chi2                  ->push_back( leadTrack->chi2()                                     );
   taus_calo_leadtrk_ndof                  ->push_back( leadTrack->ndof()                                     );
   taus_calo_leadtrk_validHits             ->push_back( leadTrack->numberOfValidHits()                        );
   taus_calo_leadtrk_lostHits              ->push_back( leadTrack->numberOfLostHits()                         );
   taus_calo_leadtrk_Signed_Sipt           ->push_back( tau_calo->leadTracksignedSipt()                       );
   taus_calo_leadtrk_HCAL3x3hitsEtSum      ->push_back( tau_calo->leadTrackHCAL3x3hitsEtSum()                 );
   taus_calo_leadtrk_HCAL3x3hottesthitDEta ->push_back( tau_calo->leadTrackHCAL3x3hottesthitDEta()            );
   taus_calo_signaltrksInvariantMass       ->push_back( tau_calo->signalTracksInvariantMass()                 );
   taus_calo_isolationtrksPtSum            ->push_back( tau_calo->isolationTracksPtSum()                      );
   taus_calo_isolationECALhitsEtSum        ->push_back( tau_calo->isolationECALhitsEtSum()                    );
   taus_calo_maximumHCALhitEt              ->push_back( tau_calo->maximumHCALhitEt()                          );
 
    
 }
  

 iEvent.put(taus_calo_p4                                ,"tauscalop4"                        );  
 iEvent.put(taus_calo_sig_ntrks                         ,"tauscalosigntrks"                  ); 
 iEvent.put(taus_calo_iso_ntrks                         ,"tauscaloisontrks"                  ); 
 iEvent.put(taus_calo_leadtrk_p4                        ,"tauscaloleadtrkp4"                 ); 
 iEvent.put(taus_calo_charge                            ,"tauscalocharge"                    ); 
 iEvent.put(taus_calo_leadtrk_d0                        ,"tauscaloleadtrkd0"                 );
 iEvent.put(taus_calo_leadtrk_z0                        ,"tauscaloleadtrkz0"                 );
 iEvent.put(taus_calo_leadtrk_chi2                      ,"tauscaloleadtrkchi2"               );
 iEvent.put(taus_calo_leadtrk_ndof                      ,"tauscaloleadtrkndof"               );
 iEvent.put(taus_calo_leadtrk_validHits                 ,"tauscaloleadtrkvalidHits"          );
 iEvent.put(taus_calo_leadtrk_lostHits                  ,"tauscaloleadtrklostHits"           );
 iEvent.put(taus_calo_leadtrk_Signed_Sipt               ,"tauscaloleadtrkSignedSipt"              );
 iEvent.put(taus_calo_leadtrk_HCAL3x3hitsEtSum          ,"tauscaloleadtrkHCAL3x3hitsEtSum"        );
 iEvent.put(taus_calo_leadtrk_HCAL3x3hottesthitDEta     ,"tauscaloleadtrkHCAL3x3hottesthitDEta"   );
 iEvent.put(taus_calo_signaltrksInvariantMass           ,"tauscalosignaltrksInvariantMass"        );
 iEvent.put(taus_calo_isolationtrksPtSum                ,"tauscaloisolationtrksPtSum"             );
 iEvent.put(taus_calo_isolationECALhitsEtSum            ,"tauscaloisolationECALhitsEtSum"         );
 iEvent.put(taus_calo_maximumHCALhitEt                  ,"tauscalomaximumHCALhitEt"               );
 
 
 
 

 
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloTauMaker);





  
