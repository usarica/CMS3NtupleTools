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
// $Id: CaloTauMaker.cc,v 1.9 2010/04/25 13:56:46 kalavase Exp $
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
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TauReco/interface/CaloTauFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoTauTag/RecoTau/interface/CaloRecoTauAlgorithm.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "CMS2/NtupleMaker/interface/CaloTauMaker.h"
#include "CMS2/NtupleMaker/interface/CommonUtils.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

CaloTauMaker::CaloTauMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_            = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<vector<vector <int> > >  (branchprefix+"isotrkidx"                     ).setBranchAlias(aliasprefix_+"_isotrk_idx"              );
  produces<vector<vector <int> > >  (branchprefix+"sigtrkidx"                     ).setBranchAlias(aliasprefix_+"_sigtrk_idx"              );
  produces<vector<int> >            (branchprefix+"leadtrkidx"                    ).setBranchAlias(aliasprefix_+"_leadtrk_idx"             );
  produces<vector<LorentzVector> >  (branchprefix+"p4"                            ).setBranchAlias(aliasprefix_+"_p4"                      );

  produces<vector<int> >            (branchprefix+"charge"                        ).setBranchAlias(aliasprefix_+"_charge"                  );
  produces<vector<float> >          (branchprefix+"leadtrkvalidHits"              ).setBranchAlias(aliasprefix_+"_leadtrk_validHits"       ); 
  produces<vector<float> >          (branchprefix+"leadtrklostHits"               ).setBranchAlias(aliasprefix_+"_leadtrk_lostHits"        ); 
  produces<vector<float> >          (branchprefix+"leadtrkSignedSipt"             ).setBranchAlias(aliasprefix_+"_leadtrk_Signed_Sipt"     );  
  produces<vector<float> >          (branchprefix+"leadtrkHCAL3x3hitsEtSum"       ).setBranchAlias(aliasprefix_+"_leadtrk_HCAL3x3hitsEtSum"); 
  produces<vector<float> >          (branchprefix+"leadtrkHCAL3x3hottesthitDEta"  ).setBranchAlias(aliasprefix_+"_leadtrk_HCAL3x3hottesthitDEta"); 
  
  produces<vector<float> >          (branchprefix+"signaltrksInvariantMass"        ).setBranchAlias(aliasprefix_+"_signaltrksInvariantMass" ); 
  produces<vector<float> >          (branchprefix+"isolationtrksPtSum"             ).setBranchAlias(aliasprefix_+"_isolationtrksPtSum"      ); 
  produces<vector<float> >          (branchprefix+"isolationECALhitsEtSum"         ).setBranchAlias(aliasprefix_+"_isolationECALhitsEtSum"  ); 
  produces<vector<float> >          (branchprefix+"maximumHCALhitEt"               ).setBranchAlias(aliasprefix_+"_maximumHCALhitEt"        ); 
  
  //tau preID
  produces<vector<int> >            (branchprefix+"tightId"                        ).setBranchAlias(aliasprefix_+"_tightId"                 );

  calotausInputTag    = iConfig.getParameter<InputTag>("calotausInputTag");
  minleadTrackPt_     = iConfig.getParameter<double>  ("minleadTrackPt");
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
  auto_ptr<vector<LorentzVector> > taus_calo_p4                               (new vector<LorentzVector> ) ;
  
  auto_ptr<vector<int> >           taus_calo_leadtrk_idx                      (new vector<int>           ) ;
  auto_ptr<vector<int> >           taus_calo_charge                           (new vector<int>) ;
  auto_ptr<vector<int> >           taus_calo_tightId                          (new vector<int>           ) ; 
  auto_ptr<vector<float> >         taus_calo_leadtrk_validHits                (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_lostHits                 (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_Signed_Sipt              (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_HCAL3x3hitsEtSum         (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_leadtrk_HCAL3x3hottesthitDEta    (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_signaltrksInvariantMass          (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_isolationtrksPtSum               (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_isolationECALhitsEtSum           (new vector<float>         ) ;
  auto_ptr<vector<float> >         taus_calo_maximumHCALhitEt                 (new vector<float>         ) ;
  
  
  
  Handle<View<reco::CaloTau> > taus_calo_h;
  iEvent.getByLabel(calotausInputTag, taus_calo_h);

  size_t tausIndex = 0;
  for(View<reco::CaloTau>::const_iterator tau_calo = taus_calo_h->begin();
      tau_calo != taus_calo_h->end(); tau_calo++, tausIndex++) {
   if(tau_calo->leadTrack().isNull())  continue;
   if(tau_calo->leadTrack()->pt()<minleadTrackPt_) continue;
   
   const edm::RefToBase<reco::CaloTau> calotauRef =  taus_calo_h->refAt(tausIndex);
   taus_calo_tightId                         ->push_back(identify(calotauRef, iSetup));
   
   const TrackRef leadTrack = tau_calo->leadTrack()  ;
   vector<int>  SigTrk_idx;
   vector<int>  IsoTrk_idx;
   const TrackRefVector& isolationTracks = tau_calo->isolationTracks();
   if( isolationTracks.size()>0){
     for(size_t iTrack = 0; iTrack < isolationTracks.size(); ++iTrack) {
       TrackRef isoTrk = isolationTracks[iTrack];
       IsoTrk_idx.push_back( static_cast<int>(isoTrk.key())  );
     }
   }

   const TrackRefVector& signalTracks = tau_calo->signalTracks();
   if( signalTracks.size()>0){
     for(size_t iTrack = 0; iTrack < signalTracks.size(); ++iTrack)
       {
	 TrackRef sigTrk = signalTracks[iTrack];
	 SigTrk_idx.push_back( static_cast<int>(sigTrk.key())  );
       }
   }

 
   taus_calo_isotrk_idx                    ->push_back( IsoTrk_idx                                              );
   taus_calo_sigtrk_idx                    ->push_back( SigTrk_idx                                              );
   taus_calo_p4                            ->push_back( LorentzVector( tau_calo->p4() )                         );
   taus_calo_charge                        ->push_back( !isfinite(tau_calo->charge()) ?
							 -9999 : tau_calo->charge());
   taus_calo_leadtrk_validHits             ->push_back( !isfinite(leadTrack->numberOfValidHits()) ? 
							-9999  :  leadTrack->numberOfValidHits()          );
							
   taus_calo_leadtrk_lostHits              ->push_back( !isfinite(leadTrack->numberOfLostHits())  ?
							-9999  : leadTrack->numberOfLostHits()            );
   taus_calo_leadtrk_Signed_Sipt           ->push_back( !isfinite(tau_calo->leadTracksignedSipt()) ? 
							-9999 : tau_calo->leadTracksignedSipt()           );
   taus_calo_leadtrk_HCAL3x3hitsEtSum      ->push_back( !isfinite(tau_calo->leadTrackHCAL3x3hitsEtSum()) ? 
							-9999 : tau_calo->leadTrackHCAL3x3hitsEtSum()     );
   taus_calo_leadtrk_HCAL3x3hottesthitDEta ->push_back( !isfinite(tau_calo->leadTrackHCAL3x3hottesthitDEta()) ?
							-9999 : tau_calo->leadTrackHCAL3x3hottesthitDEta());
   taus_calo_signaltrksInvariantMass       ->push_back( !isfinite(tau_calo->signalTracksInvariantMass()) ? 
							-9999 : tau_calo->signalTracksInvariantMass()     );
   taus_calo_isolationtrksPtSum            ->push_back( !isfinite(tau_calo->isolationTracksPtSum()) ?
							 -9999: tau_calo->isolationTracksPtSum()          );
   taus_calo_isolationECALhitsEtSum        ->push_back( !isfinite(tau_calo->isolationECALhitsEtSum()) ? 
							-9999 : tau_calo->isolationECALhitsEtSum()        );
   taus_calo_maximumHCALhitEt              ->push_back( !isfinite(tau_calo->maximumHCALhitEt()) ? 
							-9999 : tau_calo->maximumHCALhitEt()              );
   
    
 }
  
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

 iEvent.put(taus_calo_sigtrk_idx                        ,branchprefix+"sigtrkidx"                      ); 
 iEvent.put(taus_calo_isotrk_idx                        ,branchprefix+"isotrkidx"                      ); 
 iEvent.put(taus_calo_leadtrk_idx                       ,branchprefix+"leadtrkidx"                     ); 
 iEvent.put(taus_calo_p4                                ,branchprefix+"p4"                             );  
 iEvent.put(taus_calo_charge                            ,branchprefix+"charge"                         ); 
 iEvent.put(taus_calo_leadtrk_validHits                 ,branchprefix+"leadtrkvalidHits"               );
 iEvent.put(taus_calo_leadtrk_lostHits                  ,branchprefix+"leadtrklostHits"                );
 iEvent.put(taus_calo_leadtrk_Signed_Sipt               ,branchprefix+"leadtrkSignedSipt"              );
 iEvent.put(taus_calo_leadtrk_HCAL3x3hitsEtSum          ,branchprefix+"leadtrkHCAL3x3hitsEtSum"        );
 iEvent.put(taus_calo_leadtrk_HCAL3x3hottesthitDEta     ,branchprefix+"leadtrkHCAL3x3hottesthitDEta"   );
 iEvent.put(taus_calo_signaltrksInvariantMass           ,branchprefix+"signaltrksInvariantMass"        );
 iEvent.put(taus_calo_isolationtrksPtSum                ,branchprefix+"isolationtrksPtSum"             );
 iEvent.put(taus_calo_isolationECALhitsEtSum            ,branchprefix+"isolationECALhitsEtSum"         );
 iEvent.put(taus_calo_maximumHCALhitEt                  ,branchprefix+"maximumHCALhitEt"               );
 
 iEvent.put(taus_calo_tightId                           ,branchprefix+"tightId"                        );
 
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





  
