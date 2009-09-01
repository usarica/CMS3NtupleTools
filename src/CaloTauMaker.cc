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
// $Id: CaloTauMaker.cc,v 1.3 2009/09/01 07:58:43 kalavase Exp $
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
#include "FWCore/Framework/interface/EventSetup.h"

#include "CMS2/NtupleMaker/interface/CaloTauMaker.h"

#include "DataFormats/TauReco/interface/CaloTauFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoTauTag/RecoTau/interface/CaloRecoTauAlgorithm.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

CaloTauMaker::CaloTauMaker(const edm::ParameterSet& iConfig) {


  produces<vector<vector <LorentzVector> > >  ("tauscaloisotrkp4"            ).setBranchAlias("taus_calo_isotrk_p4"               );
  produces<vector<vector <LorentzVector> > >  ("tauscalosigtrkp4"            ).setBranchAlias("taus_calo_sigtrk_p4"               );
  produces<vector<LorentzVector> >  ("tauscalop4"                            ).setBranchAlias("taus_calo_p4"                      );
 //  produces<vector<int> >            ("tauscalosigntrks"                      ).setBranchAlias("taus_calo_sig_ntrks"               );
//   produces<vector<int> >            ("tauscaloisontrks"                      ).setBranchAlias("taus_calo_iso_ntrks"               );
 
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
  
  //tau preID
  produces<vector<int> >            ("tauscalotightId"                        ).setBranchAlias("taus_calo_tightId"                 );


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
  
  auto_ptr<vector<vector<LorentzVector> > > taus_calo_sigtrk_p4       (new vector<vector<LorentzVector> >) ;
  auto_ptr<vector<vector<LorentzVector> > > taus_calo_isotrk_p4       (new vector<vector<LorentzVector> >) ;
  auto_ptr<vector<LorentzVector> > taus_calo_p4                              (new vector<LorentzVector>) ;
// auto_ptr<vector<int> >           taus_calo_sig_ntrks                       (new vector<int>) ;
// auto_ptr<vector<int> >           taus_calo_iso_ntrks                       (new vector<int>) ;
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
  auto_ptr<vector<int> >           taus_calo_tightId                         (new vector<int>) ; 
  
  
  Handle<View<reco::CaloTau> > taus_calo_h;
  iEvent.getByLabel(calotausInputTag, taus_calo_h);

  
 
  size_t tausIndex = 0;
 for(View<reco::CaloTau>::const_iterator tau_calo = taus_calo_h->begin();
      tau_calo != taus_calo_h->end(); tau_calo++, tausIndex++) {
   if(tau_calo->leadTrack().isNull())  continue;
   if(tau_calo->leadTrack()->pt()<5.0) continue;
   //printf("%s  \n", "CaloTau  ");
  
   const edm::RefToBase<reco::CaloTau> calotauRef =  taus_calo_h->refAt(tausIndex);
   bool tight_id = identify(calotauRef, iSetup) ;
   if(tight_id) taus_calo_tightId                         ->push_back(1);
   else         taus_calo_tightId                         ->push_back(0);
 
   const TrackRef leadTrack = tau_calo->leadTrack()  ;
   vector<LorentzVector>  SigTrk_p4;
   vector<LorentzVector>  IsoTrk_p4;
   const TrackRefVector& isolationTracks = tau_calo->isolationTracks();
   if( isolationTracks.size()>0){
     for(size_t iTrack = 0; iTrack < isolationTracks.size(); ++iTrack)
       {
	 TrackRef isoTrk = isolationTracks[iTrack];
	 IsoTrk_p4.push_back( LorentzVector(isoTrk.get()->px(),isoTrk.get()->py(),isoTrk.get()->pz(),isoTrk.get()->p()) );
       }
   }
   else  IsoTrk_p4.push_back(LorentzVector(0, 0, 0, 0));

   const TrackRefVector& signalTracks = tau_calo->signalTracks();
   if( signalTracks.size()>0){
     for(size_t iTrack = 0; iTrack < signalTracks.size(); ++iTrack)
       {
	 TrackRef sigTrk = signalTracks[iTrack];
	 SigTrk_p4.push_back( LorentzVector(sigTrk.get()->px(),sigTrk.get()->py(),sigTrk.get()->pz(), sigTrk.get()->p()) );
       }
   }
   else  SigTrk_p4.push_back(LorentzVector(0, 0, 0, 0));

   taus_calo_isotrk_p4                     ->push_back( IsoTrk_p4                                             );
   taus_calo_sigtrk_p4                     ->push_back( SigTrk_p4                                             );
   taus_calo_p4                            ->push_back( tau_calo->p4()                                        );
//    taus_calo_sig_ntrks                     ->push_back( tau_calo->signalTracks().size()                       );
//    taus_calo_iso_ntrks                     ->push_back( tau_calo->isolationTracks().size()                    );
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
  
 iEvent.put(taus_calo_sigtrk_p4                         ,"tauscalosigtrkp4"                  ); 
 iEvent.put(taus_calo_isotrk_p4                         ,"tauscaloisotrkp4"                  ); 
 iEvent.put(taus_calo_p4                                ,"tauscalop4"                        );  
 // iEvent.put(taus_calo_sig_ntrks                         ,"tauscalosigntrks"                  ); 
//  iEvent.put(taus_calo_iso_ntrks                         ,"tauscaloisontrks"                  ); 
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
 
 iEvent.put(taus_calo_tightId                           ,"tauscalotightId"                        );
 
 }

//-------------------------------------------------------------------------------------------------
//Tau ID
//-------------------------------------------------------------------------------------------------
bool CaloTauMaker::identify(const edm::RefToBase<reco::CaloTau> &tau_calo, const edm::EventSetup& iSetup) {
 
  //isolation
  
  const TrackRefVector& isolationTracks = tau_calo->isolationTracks();
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
  //anti-electron
  bool decision = false;
  bool pass_crack = true;
  bool pass_hcaletsum = true;
  /*
  ESHandle<MagneticField> theMagneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagneticField);
  math::XYZPoint thepropagleadTrackECALSurfContactPoint=TauTagTools::propagTrackECALSurfContactPoint(theMagneticField.product(),tau_calo->leadTrack());
  if(thepropagleadTrackECALSurfContactPoint.R()==0. ||
     fabs(thepropagleadTrackECALSurfContactPoint.eta())<ECALBounds::crack_absEtaIntervalA().second || 
     (fabs(thepropagleadTrackECALSurfContactPoint.eta())>ECALBounds::crack_absEtaIntervalB().first && fabs(thepropagleadTrackECALSurfContactPoint.eta())<ECALBounds::crack_absEtaIntervalB().second) ||
     (fabs(thepropagleadTrackECALSurfContactPoint.eta())>ECALBounds::crack_absEtaIntervalC().first && fabs(thepropagleadTrackECALSurfContactPoint.eta())<ECALBounds::crack_absEtaIntervalC().second) ||
     (fabs(thepropagleadTrackECALSurfContactPoint.eta())>ECALBounds::crack_absEtaIntervalD().first && fabs(thepropagleadTrackECALSurfContactPoint.eta())<ECALBounds::crack_absEtaIntervalD().second) ||
     (fabs(thepropagleadTrackECALSurfContactPoint.eta())>ECALBounds::crack_absEtaIntervalE().first && fabs(thepropagleadTrackECALSurfContactPoint.eta())<ECALBounds::crack_absEtaIntervalE().second)
     ){
    pass_crack =  false;
  }  
  */
  float  leadTrack_HCAL3x3hitsEtSumOverPt_minvalue_ = 0.1;
  if (tau_calo->leadTrackHCAL3x3hitsEtSum()/tau_calo->leadTrack()->pt()<=leadTrack_HCAL3x3hitsEtSumOverPt_minvalue_){
    pass_hcaletsum = false;
   
  }
  decision = pass_crack && pass_hcaletsum;
  if(!decision) {
    return false;
  }
 
  return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloTauMaker);





  
