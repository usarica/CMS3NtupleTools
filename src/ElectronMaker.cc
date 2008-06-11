// -*- C++ -*-
//
// Package:    ElectronMaker
// Class:      ElectronMaker
// 
/**\class ElectronMaker ElectronMaker.cc CMS2/ElectronMaker/src/ElectronMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronMaker.cc,v 1.1 2008/06/11 03:52:33 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/ElectronMaker.h"

//#include "CMS2/Electrons/interface/Electrons.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackReco/interface/Track.h"


#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
//#include "CMS1/ExternalDataFormats/interface/EcalCluster.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;
//
// class decleration
//

//
// constructors and destructor
//
ElectronMaker::ElectronMaker(const edm::ParameterSet& iConfig)
{

  produces<vector<int> >          ("elsvalidHits"        ).setBranchAlias("els_validHits"        );
  produces<vector<int> >  	  ("elslostHits"         ).setBranchAlias("els_lostHits"         );
  produces<vector<int> >  	  ("elsmcid"             ).setBranchAlias("els_mcid"             );
  produces<vector<int> >  	  ("elscharge"           ).setBranchAlias("els_charge"           );
  produces<vector<int> >  	  ("elsmcmotherid"       ).setBranchAlias("els_mcmotherid"       );
  produces<vector<int> >  	  ("elsnSeed"            ).setBranchAlias("els_nSeed"            );
  produces<vector<int> >  	  ("elsclass"            ).setBranchAlias("els_class"            );
  produces<vector<int> >  	  ("elsrobustId"         ).setBranchAlias("els_robustId"         );
  produces<vector<int> >  	  ("elslooseId"          ).setBranchAlias("els_looseId"          );
  produces<vector<int> >  	  ("elstightId"          ).setBranchAlias("els_tightId"          );
  produces<vector<int> >  	  ("elspass3simpleId"    ).setBranchAlias("els_pass3simpleId"    );
  produces<vector<int> >	  ("elspass3looseId"     ).setBranchAlias("els_pass3looseId"     );
  produces<vector<int> >	  ("elspass3tightId"     ).setBranchAlias("els_pass3tightId"     );
  produces<vector<int> >	  ("elssimpleIdPlus"     ).setBranchAlias("els_simpleIdPlus"     );
  //produces<vector<int> >          ("elstqegammaTkNumIso" ).setBranchAlias("els_tq_egammaTkNumIso" );
  //produces<vector<int> >	    ("elstqgenID"          ).setBranchAlias("els_tq_genID"          );
  //produces<vector<int> >	    ("elstqgenMotherID"    ).setBranchAlias("els_tq_genMotherID"    );
  produces<vector<int> >	  ("elsclosestMuon"      ).setBranchAlias("els_closestMuon"      );
  produces<vector<float> >	  ("elsd0"               ).setBranchAlias("els_d0");
  produces<vector<float> >	  ("elsz0"               ).setBranchAlias("els_z0");
  produces<vector<float> >	  ("elsvertexphi"        ).setBranchAlias("els_vertexphi");
  produces<vector<float> >	  ("elschi2"             ).setBranchAlias("els_chi2");
  produces<vector<float> >	  ("elsndof"             ).setBranchAlias("els_ndof");
  produces<vector<float> >	  ("elsd0Err"            ).setBranchAlias("els_d0Err");
  produces<vector<float> >	  ("elsz0Err"            ).setBranchAlias("els_z0Err");
  produces<vector<float> >	  ("elsptErr"            ).setBranchAlias("els_ptErr");
  produces<vector<float> >	  ("elsetaErr"           ).setBranchAlias("els_etaErr");
  produces<vector<float> >	  ("elsphiErr"           ).setBranchAlias("els_phiErr");
  produces<vector<float> >	  ("elshOverE"           ).setBranchAlias("els_hOverE");
  produces<vector<float> >	  ("elseOverPIn"         ).setBranchAlias("els_eOverPIn");
  produces<vector<float> >	  ("elseSeedOverPOut"    ).setBranchAlias("els_eSeedOverPOut");
  produces<vector<float> >	  ("elsfBrem"            ).setBranchAlias("els_fBrem");
  produces<vector<float> >	  ("elsdEtaIn"           ).setBranchAlias("els_dEtaIn");
  produces<vector<float> >	  ("elsdEtaOut"          ).setBranchAlias("els_dEtaOut");
  produces<vector<float> >	  ("elsdPhiIn"           ).setBranchAlias("els_dPhiIn");
  produces<vector<float> >	  ("elsdPhiOut"          ).setBranchAlias("els_dPhiOut");
  produces<vector<float> >	  ("elsESc"              ).setBranchAlias("els_ESc");
  produces<vector<float> >	  ("elsEScraw"           ).setBranchAlias("els_ESc");
  produces<vector<float> >	  ("else3x3"             ).setBranchAlias("els_e3x3");
  produces<vector<float> >	  ("else5x5"             ).setBranchAlias("els_e5x5");
  produces<vector<float> >	  ("elsESeed"            ).setBranchAlias("els_ESeed");
  produces<vector<float> >	  ("elssigmaPhiPhi"      ).setBranchAlias("els_sigmaPhiPhi");
  produces<vector<float> >	  ("elssigmaEtaEta"      ).setBranchAlias("els_sigmaEtaEta");
  produces<vector<float> >	  ("elstkIso"            ).setBranchAlias("els_tkIso");
  //produces<vector<float> >	    ("elstqtrackIso"        ).setBranchAlias("els_tq_trackIso");
  //produces<vector<float> >	    ("elstqcaloIso"         ).setBranchAlias("els_tq_caloIso");
  //produces<vector<float> >	    ("elstqleptonID"        ).setBranchAlias("els_tq_leptonID");
  //produces<vector<float> >	    ("elstqelectronIDRobust").setBranchAlias("els_tq_electronIDRobust");
  //produces<vector<float> >	    ("elstqegammaTkIso"     ).setBranchAlias("els_tq_egammaTkIso");
  //produces<vector<float> >	    ("elstqegammaEcalIso"   ).setBranchAlias("els_tq_egammaEcalIso");
  //produces<vector<float> >	    ("elstqegammaHcalIso"   ).setBranchAlias("els_tq_egammaHcalIso");
  //produces<vector<float> >	    ("elstqLRComb"          ).setBranchAlias("els_tq_LRComb");
  produces<vector<LorentzVector> >  ("elsp4"              ).setBranchAlias("els_p4");
  produces<vector<LorentzVector> >  ("elstrkp4"           ).setBranchAlias("els_trk_p4");
  produces<vector<LorentzVector> >  ("elsmcp4"            ).setBranchAlias("els_mc_p4");
  produces<vector<LorentzVector> >  ("elsp4In"            ).setBranchAlias("els_p4_In");
  produces<vector<LorentzVector> >  ("elsp4Out"           ).setBranchAlias("els_p4_Out");
  //produces<vector<LorentzVector> >  ("elstqgenP4"         ).setBranchAlias("els_tq_genP4");
  //produces<vector<LorentzVector> >  ("elstqgenMotherP4"   ).setBranchAlias("els_tq_genMotherP4");
  
  


  //now do what ever other initialization is needed

}


ElectronMaker::~ElectronMaker()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ElectronMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  
  auto_ptr<vector<int> >           els_validHits            (new vector<int>) ;
  auto_ptr<vector<int> >	   els_lostHits             (new vector<int>) ;
  auto_ptr<vector<int> >	   els_mcid                 (new vector<int>) ;
  auto_ptr<vector<int> >	   els_charge               (new vector<int>) ;
  auto_ptr<vector<int> >	   els_mcmotherid           (new vector<int>) ;
  auto_ptr<vector<int> >	   els_nSeed                (new vector<int>) ;
  auto_ptr<vector<int> >	   els_class                (new vector<int>) ;
  auto_ptr<vector<int> >	   els_robustId             (new vector<int>) ;
  auto_ptr<vector<int> >	   els_looseId              (new vector<int>) ;
  auto_ptr<vector<int> >	   els_tightId              (new vector<int>) ;
  auto_ptr<vector<int> >	   els_pass3simpleId        (new vector<int>) ;
  auto_ptr<vector<int> >	   els_pass3looseId         (new vector<int>) ;
  auto_ptr<vector<int> >	   els_pass3tightId         (new vector<int>) ;
  auto_ptr<vector<int> >	   els_simpleIdPlus         (new vector<int>) ;
  //auto_ptr<vector<int> >	     els_tq_egammaTkNumIso    (new vector<int>) ;
  //auto_ptr<vector<int> >	     els_tq_genID             (new vector<int>) ;
  //auto_ptr<vector<int> >	     els_tq_genMotherID       (new vector<int>) ;
  auto_ptr<vector<int> >	   els_closestMuon          (new vector<int>) ;
  auto_ptr<vector<float> >	   els_d0                   (new vector<float>) ;
  auto_ptr<vector<float> >	   els_z0                   (new vector<float>) ;
  auto_ptr<vector<float> >	   els_vertexphi            (new vector<float>) ;
  auto_ptr<vector<float> >	   els_chi2                 (new vector<float>) ;
  auto_ptr<vector<float> >	   els_ndof                 (new vector<float>) ;
  auto_ptr<vector<float> >	   els_d0Err                (new vector<float>) ;
  auto_ptr<vector<float> >	   els_z0Err                (new vector<float>) ;
  auto_ptr<vector<float> >	   els_ptErr                (new vector<float>) ;
  auto_ptr<vector<float> >	   els_etaErr               (new vector<float>) ;
  auto_ptr<vector<float> >	   els_phiErr               (new vector<float>) ;
  auto_ptr<vector<float> >	   els_hOverE               (new vector<float>) ;
  auto_ptr<vector<float> >	   els_eOverPIn             (new vector<float>) ;
  auto_ptr<vector<float> >	   els_eSeedOverPOut        (new vector<float>) ;
  auto_ptr<vector<float> >	   els_fBrem                (new vector<float>) ;
  auto_ptr<vector<float> >	   els_dEtaIn               (new vector<float>) ;
  auto_ptr<vector<float> >	   els_dEtaOut              (new vector<float>) ;
  auto_ptr<vector<float> >	   els_dPhiIn               (new vector<float>) ;
  auto_ptr<vector<float> >	   els_dPhiOut              (new vector<float>) ;
  auto_ptr<vector<float> >	   els_ESc                  (new vector<float>) ;
  auto_ptr<vector<float> >	   els_EScraw               (new vector<float>) ;
  auto_ptr<vector<float> >	   els_e3x3                 (new vector<float>) ;
  auto_ptr<vector<float> >	   els_e5x5                 (new vector<float>) ;
  auto_ptr<vector<float> >	   els_ESeed                (new vector<float>) ;
  auto_ptr<vector<float> >	   els_sigmaPhiPhi          (new vector<float>) ;
  auto_ptr<vector<float> >	   els_sigmaEtaEta          (new vector<float>) ;
  auto_ptr<vector<float> >	   els_tkIso                (new vector<float>) ;
  //auto_ptr<vector<float> >         els_tq_trackIso          (new vector<float>) ;
  //auto_ptr<vector<float> >  	     els_tq_caloIso           (new vector<float>) ;
  //auto_ptr<vector<float> >  	     els_tq_leptonID          (new vector<float>) ;
  //auto_ptr<vector<float> >	     els_tq_electronIDRobust  (new vector<float>) ;
  //auto_ptr<vector<float> >	     els_tq_egammaTkIso       (new vector<float>) ;
  //auto_ptr<vector<float> >	     els_tq_egammaEcalIso     (new vector<float>) ;
  //auto_ptr<vector<float> >	     els_tq_egammaHcalIso     (new vector<float>) ;
  //auto_ptr<vector<float> >	     els_tq_LRComb            (new vector<float>) ;
  auto_ptr<vector<LorentzVector> > els_p4                   (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > els_trk_p4               (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > els_mc_p4                (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > els_p4In                 (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > els_p4Out                (new vector<LorentzVector>) ;
  //auto_ptr<vector<LorentzVector> > els_tq_genp4             (new vector<LorentzVector>) ;
  //auto_ptr<vector<LorentzVector> > els_tq_genMotherp4       (new vector<LorentzVector>) ;
  
								 
  vector<const PixelMatchGsfElectron*> electron_coll = getElectrons(iEvent);
  removeElectrons(&electron_coll);

  for(vector<const PixelMatchGsfElectron*>::const_iterator el_it = electron_coll.begin();
      el_it != electron_coll.end(); el_it++) {

    const PixelMatchGsfElectron *el = *el_it;
    const reco::Track *el_track = (const reco::Track*)(*el_it)->gsfTrack().get();
    

    els_validHits             ->push_back( el_track->numberOfValidHits()             );
    els_lostHits              ->push_back( el_track->numberOfLostHits()              );
    els_mcid                  ->push_back(999);
    els_charge                ->push_back( el->charge()                              );
    els_mcmotherid            ->push_back(999);
    els_nSeed                 ->push_back( el->numberOfClusters() - 1                );                             
    els_class                 ->push_back( el->classification()                      );
    els_robustId              ->push_back(999);
    els_looseId               ->push_back(999);
    els_tightId               ->push_back(999);
    els_pass3simpleId         ->push_back(999);
    els_pass3looseId          ->push_back(999);
    els_pass3tightId          ->push_back(99);
    els_simpleIdPlus          ->push_back(999);
    //els_tq_egammaTkNumIso     ->push_back();
    //els_tq_genID              ->push_back();
    //els_tq_genMotherID        ->push_back();
    els_closestMuon           ->push_back(999 );
    els_d0                    ->push_back( el_track->d0()                            );
    els_z0                    ->push_back( el_track->dz()                            );
    els_vertexphi             ->push_back( atan2( el_track->vy(), el_track->vx() )   );
    els_chi2                  ->push_back( el_track->chi2()                          );
    els_ndof                  ->push_back( el_track->ndof()                          );
    els_d0Err                 ->push_back( el_track->d0Error()                       );
    els_z0Err                 ->push_back( el_track->dzError()                       );
    
    float pt = el_track->pt();
    float p = el_track->p();
    float q = el_track->charge();
    float pz = el_track->pz();
    float err = (el_track->charge()!=0) ? sqrt(pt*pt*p*p/pow(q, 2)*(el_track->covariance(0,0))
                                           +2*pt*p/q*pz*(el_track->covariance(0,1))
                                           + pz*pz*(el_track->covariance(1,1) ) )
                                           : -999.;
    els_ptErr                 ->push_back( err                                       );
    els_etaErr                ->push_back( el_track->etaError()                      );
    els_phiErr                ->push_back( el_track->phiError()                      );
    els_hOverE                ->push_back( el->hadronicOverEm()                      );
    els_eOverPIn              ->push_back( el->eSuperClusterOverP()                  );
    els_eSeedOverPOut         ->push_back( el->eSeedClusterOverPout()                );

    float pin  = el->trackMomentumAtVtx().R();
    float pout = el->trackMomentumOut().R();

    els_fBrem                 ->push_back( (pin-pout)/pin                            );
    els_dEtaIn                ->push_back( el->deltaEtaSuperClusterTrackAtVtx()      );
    els_dEtaOut               ->push_back( el->deltaEtaSeedClusterTrackAtCalo()      );
    els_dPhiIn                ->push_back( el->deltaPhiSuperClusterTrackAtVtx()      );
    els_dPhiOut               ->push_back( el->deltaPhiSeedClusterTrackAtCalo()      );
    els_ESc                   ->push_back(999.);
    els_EScraw                ->push_back( el->superCluster()->rawEnergy()           );
    els_e3x3                  ->push_back(999.);
    els_e5x5                  ->push_back(999.);
    els_ESeed                 ->push_back( el->superCluster()->seed()->energy()      );
    els_sigmaPhiPhi           ->push_back(999.);
    els_sigmaEtaEta           ->push_back(999.);
    els_tkIso                 ->push_back(999.);
    //els_tq_trackIso           ->push_back();
    //els_tq_caloIso            ->push_back();
    //els_tq_leptonId           ->push_back();
    //els_tq_electronIDRobust   ->push_back();
    //els_tq_egammaTkIso        ->push_back();
    //els_tq_egammaEcalIso      ->push_back();
    //els_tq_egammaHcalIso      ->push_back();
    //els_tq_LRComb             ->push_back();
    els_p4                    ->push_back( el->p4()                                  );


    LorentzVector trk_p4( el_track->px(), el_track->py(), 
			  el_track->pz(), el_track->p() );
    
    els_trk_p4                ->push_back( trk_p4                                    );
    els_mc_p4                 ->push_back( LorentzVector(0,0,0,0)                    );
    math::XYZVector p3 = el->trackMomentumAtVtx();
    double mass = 0.000510998918;
    LorentzVector p4; 
    p4.SetXYZT(p3.x(), p3.y(), p3.z(), sqrt(mass*mass+p3.R()*p3.R()));
    
    els_p4In                  ->push_back( p4                                        );

    p3 = el->trackMomentumOut();
    p4.SetXYZT(p3.x(), p3.y(), p3.z(), sqrt(mass*mass+p3.R()*p3.R()));
        
    els_p4Out                 ->push_back( p4                                        );
    //auto_ptr<vector<LorentzVector> > els_tq_genp4             (new vector<LorentzVector>) ;
    //auto_ptr<vector<LorentzVector> > els_tq_genMotherp4       (new vector<LorentzVector>) ;
  
  }


  iEvent.put(els_validHits                    ,"elsvalidHits"     );
  iEvent.put(els_lostHits                     ,"elslostHits"      );
  iEvent.put(els_mcid                         ,"elsmcid"          );
  iEvent.put(els_charge                       ,"elscharge"        );
  iEvent.put(els_mcmotherid                   ,"elsmcmotherid"    );
  iEvent.put(els_nSeed                        ,"elsnSeed"         );
  iEvent.put(els_class                        ,"elsclass"         );
  iEvent.put(els_robustId                     ,"elsrobustId"      );
  iEvent.put(els_looseId                      ,"elslooseId"       );
  iEvent.put(els_tightId                      ,"elstightId"       );
  iEvent.put(els_pass3simpleId                ,"elspass3simpleId" );
  iEvent.put(els_pass3looseId                 ,"elspass3looseId"  );
  iEvent.put(els_pass3tightId                 ,"elspass3tightId"  );
  iEvent.put(els_simpleIdPlus                 ,"elssimpleIdPlus"  );
  //auto_ptr<vector<int> >	     els_tq_egammaTkNumIso    (new vector<int>);
  //auto_ptr<vector<int> >	     els_tq_genID             (new vector<int>);
  //auto_ptr<vector<int> >	     els_tq_genMotherID       (new vector<int>);
  iEvent.put(els_closestMuon                  ,"elsclosestMuon"          );
  iEvent.put(els_d0                   ,"elsd0"                   );
  iEvent.put(els_z0                   ,"elsz0"                   );
  iEvent.put(els_vertexphi            ,"elsvertexphi"            );
  iEvent.put(els_chi2                 ,"elschi2"                 );
  iEvent.put(els_ndof                 ,"elsndof"                 );
  iEvent.put(els_d0Err                ,"elsd0Err"                );
  iEvent.put(els_z0Err                ,"elsz0Err"                );
  iEvent.put(els_ptErr                ,"elsptErr"                );
  iEvent.put(els_etaErr               ,"elsetaErr"               );
  iEvent.put(els_phiErr               ,"elsphiErr"               );
  iEvent.put(els_hOverE               ,"elshOverE"               );
  iEvent.put(els_eOverPIn             ,"elseOverPIn"             );
  iEvent.put(els_eSeedOverPOut        ,"elseSeedOverPOut"        );
  iEvent.put(els_fBrem                ,"elsfBrem"                );
  iEvent.put(els_dEtaIn               ,"elsdEtaIn"               );
  iEvent.put(els_dEtaOut              ,"elsdEtaOut"              );
  iEvent.put(els_dPhiIn               ,"elsdPhiIn"               );
  iEvent.put(els_dPhiOut              ,"elsdPhiOut"              );
  iEvent.put(els_ESc                  ,"elsESc"                  );
  iEvent.put(els_EScraw               ,"elsEScraw"               );
  iEvent.put(els_e3x3                 ,"else3x3"                 );
  iEvent.put(els_e5x5                 ,"else5x5"                 );
  iEvent.put(els_ESeed                ,"elsESeed"                );
  iEvent.put(els_sigmaPhiPhi          ,"elssigmaPhiPhi"          );
  iEvent.put(els_sigmaEtaEta          ,"elssigmaEtaEta"          );
  iEvent.put(els_tkIso                ,"elstkIso"                );
  //iEvent.put(els_tq_trackIso          ,"elstq_trackIso"          );
  //iEvent.put(els_tq_caloIso           ,"elstq_caloIso"           );
  //iEvent.put(els_tq_leptonID          ,"elstq_leptonID"          );
  //iEvent.put(els_tq_electronIDRobust  ,"elstq_electronIDRobust"  );
  //iEvent.put(els_tq_egammaTkIso       ,"elstq_egammaTkIso"       );
  //iEvent.put(els_tq_egammaEcalIso     ,"elstq_egammaEcalIso"     );
  //iEvent.put(els_tq_egammaHcalIso     ,"elstq_egammaHcalIso"     );
  //iEvent.put(els_tq_LRComb            ,"elstq_LRComb"            );
  iEvent.put(els_p4                   ,"elsp4"                     );
  iEvent.put(els_trk_p4               ,"elstrkp4"                 );
  iEvent.put(els_mc_p4                ,"elsmcp4"                  );
  iEvent.put(els_p4In                 ,"elsp4In"                   );
  iEvent.put(els_p4Out                ,"elsp4Out"                  );
  //iEvent.put(els_tq_genp4             ,"elstq_genp4"             );
  //iEvent.put(els_tq_genMotherp4       ,"elstq_genMotherp4"       );

}




//--------------------------------------------------------------------------------------------
//gets the electrons from the event
//--------------------------------------------------------------------------------------------
vector<const PixelMatchGsfElectron*> ElectronMaker::getElectrons(const edm::Event& iEvent) {
  

  Handle<View<PixelMatchGsfElectron> > electron_h;
  //iEvent.getByLabel(electronType.c_str(), electron_h);
  iEvent.getByLabel("allLayer1TopElectrons", electron_h);
  
  vector<const PixelMatchGsfElectron*> collection;
  
  for(edm::View<reco::PixelMatchGsfElectron>::const_iterator electron = electron_h->begin(); electron != electron_h->end(); ++electron){
    collection.push_back(&*electron);
  }

  return collection;
  
}


//--------------------------------------------------------------------------------------------
//remove electrons that have the same SC
//--------------------------------------------------------------------------------------------
void ElectronMaker::removeElectrons(const vector<const PixelMatchGsfElectron*>* collection) {

  vector<const reco::PixelMatchGsfElectron*>* newEl = const_cast<vector<const PixelMatchGsfElectron*>* >(collection);
 vector<const reco::PixelMatchGsfElectron*> copy = *newEl;

  newEl->clear();

  vector<const PixelMatchGsfElectron*>::iterator it1, it2;

  for(it1=copy.begin(); it1!=copy.end(); ++it1) {

    bool isRemoved = false;
    for(it2=copy.begin(); it2!=copy.end(); ++it2) {
      if (it1 == it2)
        continue;
      if (((**it1).superCluster().id() == (**it2).superCluster().id()) &&
          ((**it1).superCluster().index() == (**it2).superCluster().index())) {

        float deltaEp1 = fabs((**it1).eSuperClusterOverP() - 1.);
        float deltaEp2 = fabs((**it2).eSuperClusterOverP() - 1.);
        if (deltaEp1 > deltaEp2) {
          isRemoved = true;
          break;
        }
      }
    }

    if (!isRemoved)
      newEl->push_back(*it1);

  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
ElectronMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMaker);
