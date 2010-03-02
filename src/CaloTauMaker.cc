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
// $Id: CaloTauMaker.cc,v 1.7 2010/03/02 19:36:07 fgolf Exp $
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

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

CaloTauMaker::CaloTauMaker(const edm::ParameterSet& iConfig) {


  produces<vector<vector <int> > >  ("tauscaloisotrkidx"                     ).setBranchAlias("taus_calo_isotrk_idx"              );
  produces<vector<vector <int> > >  ("tauscalosigtrkidx"                     ).setBranchAlias("taus_calo_sigtrk_idx"              );
  produces<vector<int> >            ("tauscaloleadtrkidx"                    ).setBranchAlias("taus_calo_leadtrk_idx"             );
  produces<vector<LorentzVector> >  ("tauscalop4"                            ).setBranchAlias("taus_calo_p4"                      );

  produces<vector<int> >            ("tauscalocharge"                        ).setBranchAlias("taus_calo_charge"                  );
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

void  CaloTauMaker::beginJob() {
}

void CaloTauMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void CaloTauMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<vector<int> > >  taus_calo_sigtrk_idx                       (new vector<vector<int> >  ) ;
  auto_ptr<vector<vector<int> > >  taus_calo_isotrk_idx                       (new vector<vector<int> >  ) ;
  auto_ptr<vector<int> >           taus_calo_leadtrk_idx                      (new vector<int>           ) ;
  auto_ptr<vector<LorentzVector> > taus_calo_p4                               (new vector<LorentzVector> ) ;
  auto_ptr<vector<int> >           taus_calo_charge                           (new vector<int>) ;

  auto_ptr<vector<float> >         taus_calo_leadtrk_validHits                (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_lostHits                 (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_Signed_Sipt              (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_HCAL3x3hitsEtSum         (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_HCAL3x3hottesthitDEta    (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_signaltrksInvariantMass          (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_isolationtrksPtSum               (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_isolationECALhitsEtSum           (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_maximumHCALhitEt                 (new vector<float>         ) ;
  auto_ptr<vector<int> >           taus_calo_tightId                          (new vector<int>           ) ; 
  
  
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
   vector<int>  SigTrk_idx;
   vector<int>  IsoTrk_idx;
   const TrackRefVector& isolationTracks = tau_calo->isolationTracks();
   if( isolationTracks.size()>0){
     for(size_t iTrack = 0; iTrack < isolationTracks.size(); ++iTrack)
       {
	 TrackRef isoTrk = isolationTracks[iTrack];
	 IsoTrk_idx.push_back( static_cast<int>(isoTrk.key())  );
       }
   }
   else  IsoTrk_idx.push_back(-9999);

   const TrackRefVector& signalTracks = tau_calo->signalTracks();
   if( signalTracks.size()>0){
     for(size_t iTrack = 0; iTrack < signalTracks.size(); ++iTrack)
       {
	 TrackRef sigTrk = signalTracks[iTrack];
	 SigTrk_idx.push_back( static_cast<int>(sigTrk.key())  );
       }
   }
   else  SigTrk_idx.push_back(-9999);

   
 
   if(leadTrack.isNonnull()){
     taus_calo_leadtrk_idx                    ->push_back( static_cast<int>(leadTrack.key())                  );
   }
   else {
         taus_calo_leadtrk_idx                ->push_back( -9999                                               );
      
   }
   taus_calo_isotrk_idx                     ->push_back( IsoTrk_idx                                           );
   taus_calo_sigtrk_idx                     ->push_back( SigTrk_idx                                           );
   taus_calo_p4                             ->push_back( LorentzVector( tau_calo->p4() )                      );
   taus_calo_charge                        ->push_back( tau_calo->charge()                                    );
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
  
 iEvent.put(taus_calo_sigtrk_idx                        ,"tauscalosigtrkidx"                      ); 
 iEvent.put(taus_calo_isotrk_idx                        ,"tauscaloisotrkidx"                      ); 
 iEvent.put(taus_calo_leadtrk_idx                       ,"tauscaloleadtrkidx"                     ); 
 iEvent.put(taus_calo_p4                                ,"tauscalop4"                             );  
 iEvent.put(taus_calo_charge                            ,"tauscalocharge"                         ); 
 iEvent.put(taus_calo_leadtrk_validHits                 ,"tauscaloleadtrkvalidHits"               );
 iEvent.put(taus_calo_leadtrk_lostHits                  ,"tauscaloleadtrklostHits"                );
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





  
