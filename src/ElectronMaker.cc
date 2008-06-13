//-*- C++ -*-
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
// $Id: ElectronMaker.cc,v 1.3 2008/06/13 19:15:20 kalavase Exp $
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

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
//#include "CMS2/ExternalDataFormats/interface/EcalCluster.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/VectorUtil.h"


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
  produces<unsigned int>            ("evtnels"             ).setBranchAlias("evt_nels"             ); //number of electrons in event
  produces<vector<int> >            ("elsvalidHits"        ).setBranchAlias("els_validHits"        ); //number of used hits in fit
  produces<vector<int> >  	    ("elslostHits"         ).setBranchAlias("els_lostHits"         ); //number of lost hits in fit
  produces<vector<int> >  	    ("elsmcid"             ).setBranchAlias("els_mcid"             ); //MC matched id
  produces<vector<int> >  	    ("elscharge"           ).setBranchAlias("els_charge"           ); //candidate charge
  produces<vector<int> >  	    ("elsmcmotherid"       ).setBranchAlias("els_mcmotherid"       ); //Id of MC matched mother
  produces<vector<int> >  	    ("elsnSeed"            ).setBranchAlias("els_nSeed"            ); 
  produces<vector<int> >  	    ("elsclass"            ).setBranchAlias("els_class"            );
  produces<vector<int> >  	    ("elsrobustId"         ).setBranchAlias("els_robustId"         );
  produces<vector<int> >  	    ("elslooseId"          ).setBranchAlias("els_looseId"          );
  produces<vector<int> >  	    ("elstightId"          ).setBranchAlias("els_tightId"          );
  produces<vector<int> >  	    ("elspass3simpleId"    ).setBranchAlias("els_pass3simpleId"    );
  produces<vector<int> >	    ("elspass3looseId"     ).setBranchAlias("els_pass3looseId"     );
  produces<vector<int> >	    ("elspass3tightId"     ).setBranchAlias("els_pass3tightId"     );
  produces<vector<int> >	    ("elssimpleIdPlus"     ).setBranchAlias("els_simpleIdPlus"     );
  produces<vector<float> >	    ("elsd0"               ).setBranchAlias("els_d0"               );
  produces<vector<float> >	    ("elsz0"               ).setBranchAlias("els_z0"               );
  produces<vector<float> >	    ("elsvertexphi"        ).setBranchAlias("els_vertexphi"        );
  produces<vector<float> >	    ("elschi2"             ).setBranchAlias("els_chi2"             );
  produces<vector<float> >	    ("elsndof"             ).setBranchAlias("els_ndof"             );
  produces<vector<float> >	    ("elsd0Err"            ).setBranchAlias("els_d0Err"            );
  produces<vector<float> >	    ("elsz0Err"            ).setBranchAlias("els_z0Err"            );
  produces<vector<float> >	    ("elsptErr"            ).setBranchAlias("els_ptErr"            );
  produces<vector<float> >	    ("elsetaErr"           ).setBranchAlias("els_etaErr"           );
  produces<vector<float> >	    ("elsphiErr"           ).setBranchAlias("els_phiErr"           );
  produces<vector<float> >	    ("elshOverE"           ).setBranchAlias("els_hOverE"           );
  produces<vector<float> >	    ("elseOverPIn"         ).setBranchAlias("els_eOverPIn"         );
  produces<vector<float> >	    ("elseSeedOverPOut"    ).setBranchAlias("els_eSeedOverPOut"    );
  produces<vector<float> >	    ("elsfBrem"            ).setBranchAlias("els_fBrem"            );
  produces<vector<float> >	    ("elsdEtaIn"           ).setBranchAlias("els_dEtaIn"           );
  produces<vector<float> >	    ("elsdEtaOut"          ).setBranchAlias("els_dEtaOut"          );
  produces<vector<float> >	    ("elsdPhiIn"           ).setBranchAlias("els_dPhiIn"           );
  produces<vector<float> >	    ("elsdPhiOut"          ).setBranchAlias("els_dPhiOut"          );
  produces<vector<float> >	    ("elsESc"              ).setBranchAlias("els_ESc"              );
  produces<vector<float> >	    ("elsEScraw"           ).setBranchAlias("els_ESc_raw"          );
  produces<vector<float> >	    ("else3x3"             ).setBranchAlias("els_e3x3"             );
  produces<vector<float> >	    ("else5x5"             ).setBranchAlias("els_e5x5"             );
  produces<vector<float> >	    ("elsESeed"            ).setBranchAlias("els_ESeed"            );
  produces<vector<float> >	    ("elssigmaPhiPhi"      ).setBranchAlias("els_sigmaPhiPhi"      );
  produces<vector<float> >	    ("elssigmaEtaEta"      ).setBranchAlias("els_sigmaEtaEta"      );
  produces<vector<float> >	    ("elstkIso"            ).setBranchAlias("els_tkIso"            );
  produces<vector<float> >	    ("elsdPhiInPhiOut"     ).setBranchAlias("els_dPhiInPhiOut"     );
  //produces<vector<float> >	    ("elsdPhiOutEcalHit"   ).setBranchAlias("els_dPhiOutEcalHit"   );
  produces<vector<LorentzVector> >  ("elsp4"               ).setBranchAlias("els_p4"               );
  produces<vector<LorentzVector> >  ("elstrkp4"            ).setBranchAlias("els_trk_p4"           );
  produces<vector<LorentzVector> >  ("elsmcp4"             ).setBranchAlias("els_mc_p4"            );
  produces<vector<LorentzVector> >  ("elsp4In"             ).setBranchAlias("els_p4In"            );
  produces<vector<LorentzVector> >  ("elsp4Out"            ).setBranchAlias("els_p4Out"           );
  
}


ElectronMaker::~ElectronMaker()
{
}

void  ElectronMaker::beginJob(const edm::EventSetup&)
{
}

void ElectronMaker::endJob() {
}



// ------------ method called to produce the data  ------------
void ElectronMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  

  auto_ptr<unsigned int>           evt_nels                 (new unsigned int         ) ;
  auto_ptr<vector<int> >           els_validHits            (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_lostHits             (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_mcid                 (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_charge               (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_mcmotherid           (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_nSeed                (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_class                (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_robustId             (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_looseId              (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_tightId              (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_pass3simpleId        (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_pass3looseId         (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_pass3tightId         (new vector<int>          ) ;
  auto_ptr<vector<int> >	   els_simpleIdPlus         (new vector<int>          ) ;
  auto_ptr<vector<float> >	   els_d0                   (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_z0                   (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_vertexphi            (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_chi2                 (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_ndof                 (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_d0Err                (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_z0Err                (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_ptErr                (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_etaErr               (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_phiErr               (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_hOverE               (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_eOverPIn             (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_eSeedOverPOut        (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_fBrem                (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_dEtaIn               (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_dEtaOut              (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_dPhiIn               (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_dPhiOut              (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_ESc                  (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_EScraw               (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_e3x3                 (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_e5x5                 (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_ESeed                (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_sigmaPhiPhi          (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_sigmaEtaEta          (new vector<float>        ) ;
  auto_ptr<vector<float> >	   els_tkIso                (new vector<float>        ) ;
  auto_ptr<vector<float> >         els_dPhiInPhiOut         (new vector<float>        ) ;
  //auto_ptr<vector<float> >	   els_dPhiOutEcalHit       (new vector<float>        ) ;
  auto_ptr<vector<LorentzVector> > els_p4                   (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > els_trk_p4               (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > els_mc_p4                (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > els_p4In                 (new vector<LorentzVector>) ;
  auto_ptr<vector<LorentzVector> > els_p4Out                (new vector<LorentzVector>) ;


  //Get required collections ->have to be made configurable, hard wired for now
  vector<const PixelMatchGsfElectron*> electron_coll = getElectrons(iEvent);
  removeElectrons(&electron_coll);

  Handle<View<Track> > tk_h;
  iEvent.getByLabel("ctfWithMaterialTracks", tk_h);
  const View<Track> *track_coll = tk_h.product();
  
  //for(View<Track>::const_iterator tk_it = tk_h->begin();
  //tk_it != tk_h->end(); tk_it++) {
  //track_coll->push_back(*tk_it);
  //}
  
  Handle<BasicClusterShapeAssociationCollection> barrelClShp_h;
  iEvent.getByLabel("hybridSuperClusters", "hybridShapeAssoc", barrelClShp_h);
  if(!barrelClShp_h.isValid()) {
    LogError("ElectronMaker") << "Can't get EcalBarrel Cluster Shape Collection";
  }
  const BasicClusterShapeAssociationCollection *barrelClShp = barrelClShp_h.product();
  
  Handle<BasicClusterShapeAssociationCollection> endcapClShp_h;
  iEvent.getByLabel("islandBasicClusters", "islandEndcapShapeAssoc", endcapClShp_h);
  if(!barrelClShp_h.isValid()) {
    LogError("ElectronMaker") << "Can't get Ecal EndcapCluster Shape Collection";
  }
  const BasicClusterShapeAssociationCollection *endcapClShp = endcapClShp_h.product();
  
  //fill number of electrons variable
  *evt_nels = electron_coll.size();
  
  //loop over electron collection
  for(vector<const PixelMatchGsfElectron*>::const_iterator el_it = electron_coll.begin();
      el_it != electron_coll.end(); el_it++) {

    const PixelMatchGsfElectron *el = *el_it;
    const reco::Track *el_track = (const reco::Track*)(*el_it)->gsfTrack().get();
    

    //get Cluster info
    //const EcalSuperClusterCollection *barrelExtraShapes, *endcapExtraShapes;
    float esc, e3x3, e5x5, see, spp;
    R9_25(barrelClShp, endcapClShp, el, esc, e3x3, e5x5, spp, see);
    

    
    /* fill Electron Id array
       robustId = 0, looseId = 1, tightId = 2, simpleId = 3,
       oldlooseId = 4, oldtightId = 5, simpleIdPlus = 6
    */
    int id[7];
    for(int i=0; i<7; ++i) {
      if (identify(el,barrelClShp,endcapClShp, i))
        id[i] = 1;
      else
        id[i] = 0;
    }


    //stuff for track pterr
    float pt = el_track->pt();
    float p = el_track->p();
    float q = el_track->charge();
    float pz = el_track->pz();
    float err = (el_track->charge()!=0) ? sqrt(pt*pt*p*p/pow(q, 2)*(el_track->covariance(0,0))
					       +2*pt*p/q*pz*(el_track->covariance(0,1))
					       + pz*pz*(el_track->covariance(1,1) ) ) : -999.;
    
    //variables to get fbrem
    float pin  = el->trackMomentumAtVtx().R();
    float pout = el->trackMomentumOut().R();

    
    //electron track p4
    LorentzVector trk_p4( el_track->px(), el_track->py(), 
			  el_track->pz(), el_track->p() );
    

    //inner and outer p4s
    double          mass = 0.000510998918;
    LorentzVector   p4In; 
    math::XYZVector p3In = el->trackMomentumAtVtx();
    p4In.SetXYZT(   p3In.x(), p3In.y(), p3In.z(), sqrt(mass*mass+p3In.R()*p3In.R()));
    
    LorentzVector   p4Out; 
    math::XYZVector p3Out = el->trackMomentumOut();
    p4Out.SetXYZT(  p3Out.x(), p3Out.y(), p3Out.z(), sqrt(mass*mass+p3Out.R()*p3Out.R()));

    
    float tkIso = trackRelIsolation(el->momentum(), el->vertex(), track_coll,
				    0.3,     //! dR < 0.3
				    0.01,    //! dR > 0.01
				    0.1,     //! |d0_tk| < 0.1 cm
				    999.9,   //! |el_2D - tk_2D| < 999
				    0.5,     //! |z0_el - z0_track| < 0.5
				    1.0,     //! min pt
				    7);      //! min nHits

    //variables need to fill dPhiInPhiOut
    double phi_pin = el->caloPosition().phi() -
      el->deltaPhiSuperClusterTrackAtVtx();
    double phi_pout = el->superCluster()->seed()->position().phi()
      - el->deltaPhiSeedClusterTrackAtCalo();


    //for dPhiOutEcalHit
    //  for (EcalSuperClusterCollection::const_iterator i_barrel = barrelExtraShapes->begin();
    // 	 i_barrel != barrelExtraShapes->end(); ++i_barrel) {
    //       if (i_barrel->superCluster() == aElectron->superCluster()) 
    // 	break;
    //     }

    
    //     for (EcalSuperClusterCollection::const_iterator i_endcap = endcapExtraShapes->begin();
    // 	 i_endcap != endcapExtraShapes->end(); ++i_endcap) {
    //       if (i_endcap->superCluster() == aElectron->superCluster()) 
    // 	   break;
    //     }
    
    //     float dPhiOutEcalHit = -999.;
    //     if (i_barrel != barrelExtraShapes->end()) {
    //       dPhiOutEcalHit = i_barrel->DphitoNonzeroCrystal(phi_pout);
    //     } else if (i_endcap != endcapExtraShapes->end()) {
    //       dPhiOutEcalHit = i_endcap->DphitoNonzeroCrystal(phi_pout);
    //    }

    
    //fill the vectors
    els_validHits             ->push_back( el_track->numberOfValidHits()             );
    els_lostHits              ->push_back( el_track->numberOfLostHits()              );
    els_mcid                  ->push_back(999); //placeholder 
    els_charge                ->push_back( el->charge()                              );
    els_mcmotherid            ->push_back(999); //placeholder 
    els_nSeed                 ->push_back( el->numberOfClusters() - 1                );                             
    els_class                 ->push_back( el->classification()                      );
    els_robustId              ->push_back( id[0]                                     );
    els_looseId               ->push_back( id[1]                                     );
    els_tightId               ->push_back( id[2]                                     );
    els_pass3simpleId         ->push_back( id[3]                                     );
    els_pass3looseId          ->push_back( id[4]                                     );
    els_pass3tightId          ->push_back( id[5]                                     );
    els_simpleIdPlus          ->push_back( id[6]                                     );
    els_d0                    ->push_back( el_track->d0()                            );
    els_z0                    ->push_back( el_track->dz()                            );
    els_vertexphi             ->push_back( atan2( el_track->vy(), el_track->vx() )   );
    els_chi2                  ->push_back( el_track->chi2()                          );
    els_ndof                  ->push_back( el_track->ndof()                          );
    els_d0Err                 ->push_back( el_track->d0Error()                       );
    els_z0Err                 ->push_back( el_track->dzError()                       );
    els_ptErr                 ->push_back( err                                       );
    els_etaErr                ->push_back( el_track->etaError()                      );
    els_phiErr                ->push_back( el_track->phiError()                      );
    els_hOverE                ->push_back( el->hadronicOverEm()                      );
    els_eOverPIn              ->push_back( el->eSuperClusterOverP()                  );
    els_eSeedOverPOut         ->push_back( el->eSeedClusterOverPout()                );
    els_fBrem                 ->push_back( (pin-pout)/pin                            );
    els_dEtaIn                ->push_back( el->deltaEtaSuperClusterTrackAtVtx()      );
    els_dEtaOut               ->push_back( el->deltaEtaSeedClusterTrackAtCalo()      );
    els_dPhiIn                ->push_back( el->deltaPhiSuperClusterTrackAtVtx()      );
    els_dPhiOut               ->push_back( el->deltaPhiSeedClusterTrackAtCalo()      );
    els_ESc                   ->push_back( esc                                       );
    els_EScraw                ->push_back( el->superCluster()->rawEnergy()           );
    els_e3x3                  ->push_back( e3x3                                      );
    els_e5x5                  ->push_back( e5x5                                      );
    els_ESeed                 ->push_back( el->superCluster()->seed()->energy()      );
    els_sigmaPhiPhi           ->push_back( see                                       );
    els_sigmaEtaEta           ->push_back( spp                                       );
    els_tkIso                 ->push_back( tkIso                                     );
    els_dPhiInPhiOut          ->push_back( phi_pin - phi_pout                        );
    //els_dPhiOutEcalHit        ->push_back( dPhiOutEcalHit                            );
    els_p4                    ->push_back( el->p4()                                  );
    els_trk_p4                ->push_back( trk_p4                                    );
    els_mc_p4                 ->push_back( LorentzVector(0,0,0,0)                    );
    els_p4In                  ->push_back( p4In                                      );
    els_p4Out                 ->push_back( p4Out                                     );
  
  }

  iEvent.put(evt_nels                         ,"evtnels"            );
  iEvent.put(els_validHits                    ,"elsvalidHits"       );
  iEvent.put(els_lostHits                     ,"elslostHits"        );
  iEvent.put(els_mcid                         ,"elsmcid"            );
  iEvent.put(els_charge                       ,"elscharge"          );
  iEvent.put(els_mcmotherid                   ,"elsmcmotherid"      );
  iEvent.put(els_nSeed                        ,"elsnSeed"           );
  iEvent.put(els_class                        ,"elsclass"           );
  iEvent.put(els_robustId                     ,"elsrobustId"        );
  iEvent.put(els_looseId                      ,"elslooseId"         );
  iEvent.put(els_tightId                      ,"elstightId"         );
  iEvent.put(els_pass3simpleId                ,"elspass3simpleId"   );
  iEvent.put(els_pass3looseId                 ,"elspass3looseId"    );
  iEvent.put(els_pass3tightId                 ,"elspass3tightId"    );
  iEvent.put(els_simpleIdPlus                 ,"elssimpleIdPlus"    );
  iEvent.put(els_d0                           ,"elsd0"              );
  iEvent.put(els_z0                           ,"elsz0"              );
  iEvent.put(els_vertexphi                    ,"elsvertexphi"       );
  iEvent.put(els_chi2                         ,"elschi2"            );
  iEvent.put(els_ndof                         ,"elsndof"            );
  iEvent.put(els_d0Err                        ,"elsd0Err"           );
  iEvent.put(els_z0Err                        ,"elsz0Err"           );
  iEvent.put(els_ptErr                        ,"elsptErr"           );
  iEvent.put(els_etaErr                       ,"elsetaErr"          );
  iEvent.put(els_phiErr                       ,"elsphiErr"          );
  iEvent.put(els_hOverE                       ,"elshOverE"          );
  iEvent.put(els_eOverPIn                     ,"elseOverPIn"        );
  iEvent.put(els_eSeedOverPOut                ,"elseSeedOverPOut"   );
  iEvent.put(els_fBrem                        ,"elsfBrem"           );
  iEvent.put(els_dEtaIn                       ,"elsdEtaIn"          );
  iEvent.put(els_dEtaOut                      ,"elsdEtaOut"         );
  iEvent.put(els_dPhiIn                       ,"elsdPhiIn"          );
  iEvent.put(els_dPhiOut                      ,"elsdPhiOut"         );
  iEvent.put(els_ESc                          ,"elsESc"             );
  iEvent.put(els_EScraw                       ,"elsEScraw"          );
  iEvent.put(els_e3x3                         ,"else3x3"            );
  iEvent.put(els_e5x5                         ,"else5x5"            );
  iEvent.put(els_ESeed                        ,"elsESeed"           );
  iEvent.put(els_sigmaPhiPhi                  ,"elssigmaPhiPhi"     );
  iEvent.put(els_sigmaEtaEta                  ,"elssigmaEtaEta"     );
  iEvent.put(els_tkIso                        ,"elstkIso"           );
  iEvent.put(els_dPhiInPhiOut                 ,"elsdPhiInPhiOut"    );
  //iEvent.put(els_dPhiOutEcalHit               ,"elsdPhiOutEcalHit"  );
  iEvent.put(els_p4                           ,"elsp4"              );
  iEvent.put(els_trk_p4                       ,"elstrkp4"           );
  iEvent.put(els_mc_p4                        ,"elsmcp4"            );
  iEvent.put(els_p4In                         ,"elsp4In"            );
  iEvent.put(els_p4Out                        ,"elsp4Out"           );

}




//--------------------------------------------------------------------------------------------
//gets the electrons from the event
//--------------------------------------------------------------------------------------------
vector<const PixelMatchGsfElectron*> ElectronMaker::getElectrons(const edm::Event& iEvent) {
  

  Handle<View<PixelMatchGsfElectron> > electron_h;
  //iEvent.getByLabel(electronType.c_str(), electron_h);
  iEvent.getByLabel("allLayer1TopElectrons", electron_h);
  
  vector<const PixelMatchGsfElectron*> collection;
  
  for(edm::View<reco::PixelMatchGsfElectron>::const_iterator electron = electron_h->begin(); 
      electron != electron_h->end(); ++electron){
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

//--------------------------------------------------------------------------------------------
//Function to fill SC Parameters
//--------------------------------------------------------------------------------------------
void ElectronMaker::R9_25(const reco::BasicClusterShapeAssociationCollection* barrelClShp,
			  const reco::BasicClusterShapeAssociationCollection* endcapClShp, 
			  const reco::PixelMatchGsfElectron* electron,
			  float& eMax, float& e3x3, float& e5x5, float& spp, float& see) { 

  reco::SuperClusterRef sclRef=electron->superCluster();
  reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
  DetId id = sclRef->seed()->getHitsByDetId()[0];
  if (id.subdetId() == EcalBarrel) {
    seedShpItr = barrelClShp->find(sclRef->seed());
  } else {
    seedShpItr = endcapClShp->find(sclRef->seed());
  }
  
  // Get the ClusterShapeRef corresponding to the BasicCluster
  const reco::ClusterShapeRef& seedShapeRef = seedShpItr->val;
  
  eMax = sclRef->energy();
  e3x3 = seedShapeRef->e3x3();
  e5x5 = seedShapeRef->e5x5();
  spp = sqrt(seedShapeRef->covPhiPhi());
  see = sqrt(seedShapeRef->covEtaEta());
}


//-------------------------------------------------------------------------------------------------
//Electron ID
//-------------------------------------------------------------------------------------------------
bool ElectronMaker::identify(const reco::PixelMatchGsfElectron* electron, 
			     const reco::BasicClusterShapeAssociationCollection* barrelClShp,
			     const reco::BasicClusterShapeAssociationCollection* endcapClShp, int type) {
  
  // type = 0 robust (looser and safe), 1 loose, 2 tight, 3 simple (old robust), 4 old loose (6 categories), 5 old tight (6 categories) 
  
  float dummy, e3x3, e5x5, i, l; 
  R9_25(barrelClShp, endcapClShp, electron, dummy, e3x3, e5x5, l, i);
  int eb;
  
  double eOverP = electron->eSuperClusterOverP();
  double eSeedOverPOut = electron->eSeedClusterOverPout();
  double eSeed = electron->superCluster()->seed()->energy();
  double pin  = electron->trackMomentumAtVtx().R();   
  double eSeedOverPin = eSeed/pin; 
  double pout = electron->trackMomentumOut().R(); 
  double fBrem = (pin-pout)/pin;
  double hOverE = electron->hadronicOverEm();
  double sigmaee = i;
  double sigmaeecor = i;
  double deltaPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();
  double deltaEtaIn = electron->deltaEtaSuperClusterTrackAtVtx();
  double eta = fabs(electron->p4().Eta());
  double d0 = electron->gsfTrack()->d0();

  if (eta < 1.479) 
    eb = 0;
  else {
    eb = 1; 
    sigmaeecor = sigmaee - 0.02*(fabs(eta) - 2.3);   //correct sigmaetaeta dependence on eta in endcap
  }
  
  if (type == 2) { 
    
    float cuthoe[]  = {0.05,    0.042,  0.045,  0.,  
                       0.055,   0.037,  0.050,  0.};
    float cutsee[]  = {0.0125,  0.011,  0.01,   0.,
                       0.0265,  0.0252, 0.026,  0.};
    float cutdphi[] = {0.032,   0.016,  0.0525, 0.09, 
                       0.025,   0.035,  0.065,  0.092};
    float cutdeta[] = {0.0055,  0.0030, 0.0065, 0.,
                       0.0060,  0.0055, 0.0075, 0.}; 
    float cuteop2[] = {0.24,    0.94,   0.11,   0.,   
                       0.32,    0.83,   0.,     0.}; 
    
    int cat = classify(electron);
        
    // TIGHT Selection
    if ((eOverP < 0.8) && (fBrem < 0.2)) 
      return false;
    
    if (eOverP < 0.9*(1-fBrem))
      return false;
    
    if (hOverE > cuthoe[cat+eb*4]) 
      return false;    
    
    if (sigmaeecor > cutsee[cat+4*eb]) 
      return false;    
    
    if (eOverP < 1.5) {
      if (fabs(deltaPhiIn) > cutdphi[cat+4*eb]) 
        return false;    
    } else {
      if (fabs(deltaPhiIn) > cutdphi[3+4*eb])
        return false;
    }
    
    if (fabs(deltaEtaIn) > cutdeta[cat+4*eb]) 
      return false;    
    
    if (eSeedOverPin < cuteop2[cat+4*eb]) 
      return false;  
    
    return true; 
  }

  if (type == 1) { 
    
    float cutsee[]  = {0.014,  0.012,  0.0115, 0.,
                       0.0275, 0.0265, 0.0265, 0.};  // first row barrel, second endcap
    float cutdeta[] = {0.009,  0.0045, 0.0085, 0.,
                       0.0105, 0.0068, 0.010,  0.};
    float cuteop2[] = {0.11,   0.91,    0.11,  0.,
                       0.0,    0.85,    0.0,   0.};
    float cutdphi[] = {0.05,   0.025,  0.053,  0.09,
                       0.07,   0.03,   0.092,  0.092};
    float cuthoe[]  = {0.115,  0.10,   0.055,  0.,
                       0.145,  0.12,   0.15,   0.};
    
    int cat = classify(electron);
    
    // const reco::ClusterShapeRef& shapeRef = getClusterShape(electron, e);
    
    // LOOSE Selection
    if ((eOverP < 0.8) && (fBrem < 0.2)) 
      return false;
    
    if (hOverE > cuthoe[cat+eb*4]) 
      return false;    
    
    if (sigmaeecor > cutsee[cat+4*eb]) 
      return false;    
    
    if (eOverP < 1.5) {
      if (fabs(deltaPhiIn) > cutdphi[cat+4*eb]) 
	return false;    
    } else {
      if (fabs(deltaPhiIn) > cutdphi[3+4*eb])
	return false;
    }
    
    if (fabs(deltaEtaIn) > cutdeta[cat+4*eb]) 
      return false;    
    
    if (eSeedOverPin < cuteop2[cat+4*eb]) 
      return false;  
    
    return true; 
  }

  if (type == 0) {
    
    double cut[] = {0.115, 0.0140, 0.090, 0.0090,
                    0.150, 0.0275, 0.092, 0.0105};
    
    if (hOverE > cut[0+eb*4]) 
      return false;    
    
    if (sigmaeecor > cut[1+eb*4]) 
      return false;    
    
    if (fabs(deltaPhiIn) > cut[2+eb*4]) 
      return false;    
    
    if (fabs(deltaEtaIn) > cut[3+eb*4]) 
      return false;    
    
    return true;
  }
  
  if (type == 3) {
    double dRIn = sqrt(pow(deltaEtaIn,2)+pow(deltaPhiIn,2));
    if ((eta > 1.479 || sigmaee < 0.012) && sigmaee < 0.03 && 
        fabs(deltaEtaIn) < 0.012 && hOverE < 0.1 && eOverP > 0.9 && dRIn < 0.1)
      return true;
    
    return false;
  }
  
  if (type == 4) {
    
    int cat_old = classify_old(electron);
    
    float cutsee[] =  {0.014,  0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 
                       0.031,  0.028,  0.028,  0.028,  0.028,  0.028,  0.028,  0.028}; // first row barrel, second endcap
    float cutdeta[] = {0.009,  0.0085, 0.0085, 0.0035, 0.0085, 0.0035, 0.0085, 0.0085,
                       0.009,  0.009,  0.009,  0.009,  0.009,  0.009,  0.009,  0.009};
    float cutdphi[] = {0.06,   0.06,   0.10,   0.02,   0.06,   0.06,   0.06,   0.06,
                       0.06,   0.10,   0.10,   0.02,   0.10,   0.10,   0.10,   0.10};
    float cuthoe[] =  {0.125,  0.055,  0.055,  0.10,   0.055,  0.055,  0.055,  0.055,
                       0.125,  0.10,   0.10,   0.10,   0.10,   0.10,   0.10,   0.10};
    float cuteop2[] = {0.55,   0.88,   0.88,   0.88,   0.88,   0.88,   0.88,   0.88,
                       0.85,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85};
    
    if(sigmaeecor < cutsee[cat_old+eb*8]) {
      if(fabs(deltaEtaIn) < cutdeta[cat_old+eb*8]) {
        if(fabs(deltaPhiIn) < cutdphi[cat_old+eb*8]) {
          if(eSeedOverPOut > cuteop2[cat_old+eb*8]) {
            if(hOverE < cuthoe[cat_old+eb*8]) {
              return true;
            }
          }
        }
      }
    }
  
    return false;
  }

  if (type == 5) {
    
    int cat_old = classify_old(electron);
    
    float cutsee[] =  {0.014,  0.012,  0.012,  0.012,  0.012,  0.012,  0.012,  0.012, 
                       0.027,  0.026,  0.026,  0.026,  0.026,  0.026,  0.026,  0.026}; // first row barrel, second endcap
    float cutdeta[] = {0.009,  0.0055, 0.0055, 0.004,  0.0055, 0.004,  0.0055, 0.0055,
                       0.008,  0.008,  0.008,  0.008,  0.008,  0.008,  0.008,  0.008};
    float cutdphi[] = {0.04,   0.04,   0.07,   0.02,   0.04,   0.04,   0.04,   0.04,
                       0.055,  0.085,  0.085,  0.02,   0.085,  0.085,  0.085,  0.085};
    float cuthoe[] =  {0.12,   0.055,  0.055,  0.095,  0.055,  0.055,  0.055,  0.055,
                       0.05,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03};
    float cuteop2[] = {0.55,   0.91,   0.91,   0.91,   0.91,   0.91,   0.91,   0.91,
                       0.85,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85};
    
    if(sigmaeecor < cutsee[cat_old+eb*8]) {
      if(fabs(deltaEtaIn) < cutdeta[cat_old+eb*8]) {
        if(fabs(deltaPhiIn) < cutdphi[cat_old+eb*8]) {
          if(eSeedOverPOut > cuteop2[cat_old+eb*8]) {
            if(hOverE < cuthoe[cat_old+eb*8]) {
              return true;
            }
          }
        }
      }
    }
    
    return false;
  }

  //This is for the simple selection Plus
  //Added by mnorman using types suggested by A. Yagil
  if (type == 6) {
    double dRIn = sqrt(pow(deltaEtaIn,2)+pow(deltaPhiIn,2));
    if ((eta > 1.479 || sigmaee < 0.012) && sigmaee < 0.03 && 
        fabs(deltaEtaIn) < 0.012 && hOverE < 0.05 && eOverP > 0.9 && dRIn < 0.1 && d0 < 0.02)
      return true;
    
    return false;
  }
  
  return false;
}

//----------------------------------------------------------------------------
//Old Electron Id classification function
//----------------------------------------------------------------------------
int ElectronMaker::classify_old(const reco::PixelMatchGsfElectron* electron) {
   
  double eOverP = electron->eSuperClusterOverP();
  double pin  = electron->trackMomentumAtVtx().R(); 
  double pout = electron->trackMomentumOut().R(); 
  double fBrem = (pin-pout)/pin;
  
  int cat=6;
  if(eOverP > 0.9*(1-fBrem)) {
    if(fBrem > 0.06) {
      if(eOverP > 1.5) 
        cat=2;
      else if(eOverP > 1.2) 
        cat=1;
      else if(eOverP > 0.8) 
        cat=0;
      else 
        cat=4;
    }
    else if(eOverP > 1.5) 
      cat=5;
    else 
      cat=3;
  } 
  
  return cat;
}

//----------------------------------------------------------------------------
//Electron Id classification function
//----------------------------------------------------------------------------
int ElectronMaker::classify(const reco::PixelMatchGsfElectron* electron) {
  
  double eta = electron->p4().Eta();
  double eOverP = electron->eSuperClusterOverP();
  double pin  = electron->trackMomentumAtVtx().R(); 
  double pout = electron->trackMomentumOut().R(); 
  double fBrem = (pin-pout)/pin;
  
  int cat;
  
  if((fabs(eta)<1.479 && fBrem<0.06) || (fabs(eta)>1.479 && fBrem<0.1)) 
    cat=1;
  else if (eOverP < 1.2 && eOverP > 0.8) 
    cat=0;
  else 
    cat=2;
  
  return cat;
}


//---------------------------------------------------------------------------
//Track Isolation
//---------------------------------------------------------------------------
double ElectronMaker::trackRelIsolation(const math::XYZVector momentum,
					const math::XYZPoint vertex,
					const  edm::View<reco::Track>* tracks,
					double dRConeMax, double dRConeMin,
					double tkVtxDMax,
					double vtxDiffDMax, double vtxDiffZMax, double ptMin, unsigned int nHits)
{
  double isoResult = -10.;
  if ( tracks == 0 ) {
    std::cout << "Configuration Error: track collection is not set!" <<std::endl;
    return isoResult;
  }

  double sumPt = 0;

  View<Track>::const_iterator iTk = tracks->begin();
  for (; iTk != tracks->end(); ++iTk){
    double dR = ROOT::Math::VectorUtil::DeltaR(momentum, iTk->momentum());
    //exclude tks in veto cone (set it to small number to
    //exclude this track
    double dZ = fabs(vertex.z() - iTk->vz());
    double d0 = sqrt(iTk->vertex().perp2());
    double dD0 = sqrt((iTk->vertex() - vertex).perp2());
    if (dR < dRConeMin) continue;
    if ( dR < dRConeMax
	 && dZ < vtxDiffZMax
	 && d0 < tkVtxDMax
	 && dD0 < vtxDiffDMax
	 && iTk->pt() >= ptMin
	 && iTk->found() > nHits){
      sumPt += iTk->pt();
    }
  }

  isoResult = sumPt;
  return isoResult;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMaker);
