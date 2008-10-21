// -*- C++ -*-
//
// Package:    MuonMaker
// Class:      MuonMaker
// 
/**\class MuonMaker MuonMaker.cc CMS2/MuonMaker/src/MuonMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonMaker.cc,v 1.12 2008/10/21 16:27:55 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/MuonMaker.h"
#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
//#include "RecoMuon/MuonIdentification/interface/IdGlobalFunctions.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Math/VectorUtil.h"


typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;
using namespace std;
using namespace reco;

MuonMaker::MuonMaker(const edm::ParameterSet& iConfig)
{
  // mu track quantities
  produces<vector<int> >	     ("mustype"		).setBranchAlias("mus_type"          );	// type
  produces<vector<int> >	     ("musgoodmask"	).setBranchAlias("mus_goodmask"      ); // good mask
  produces<vector<LorentzVector> >   ("musp4"		).setBranchAlias("mus_p4"            ); // candidate p4						
  produces<vector<LorentzVector> >   ("mustrkp4"	).setBranchAlias("mus_trk_p4"        ); // track p4						
  produces<vector<float> >	     ("musd0"		).setBranchAlias("mus_d0"            ); // impact parameter at the point of closest approach	
  produces<vector<float> >	     ("musz0"		).setBranchAlias("mus_z0"            ); // z position of the point of closest approach		
  produces<vector<float> >	     ("musd0corr"	).setBranchAlias("mus_d0corr"        ); // corrected impact parameter at the point of closest approach	
  produces<vector<float> >	     ("musz0corr"	).setBranchAlias("mus_z0corr"        ); // corrected z position of the point of closest approach		
  produces<vector<float> >	     ("musvertexphi"	).setBranchAlias("mus_vertexphi"     ); // phi angle of the point of closest approach		
  produces<vector<float> >	     ("muschi2"		).setBranchAlias("mus_chi2"          ); // chi2 of the silicon tracker fit			
  produces<vector<float> >	     ("musndof"		).setBranchAlias("mus_ndof"          ); // number of degrees of freedom of the fit		
  produces<vector<int> >	     ("musvalidHits"	).setBranchAlias("mus_validHits"     ); // number of used hits in the fit			
  produces<vector<int> >	     ("muslostHits"	).setBranchAlias("mus_lostHits"      ); // number of lost hits in the fit			
  produces<vector<float> >	     ("musd0Err"	).setBranchAlias("mus_d0Err"         ); // error on the impact parameter			
  produces<vector<float> >	     ("musz0Err"	).setBranchAlias("mus_z0Err"         ); // error on z position of the point of closest approach	
  produces<vector<float> >	     ("musptErr"	).setBranchAlias("mus_ptErr"         ); // track Pt error					
  produces<vector<float> >	     ("musetaErr"	).setBranchAlias("mus_etaErr"        ); // track eta error					
  produces<vector<float> >	     ("musphiErr"	).setBranchAlias("mus_phiErr"        ); // track phi error					
  produces<vector<int> >	     ("muscharge"	).setBranchAlias("mus_charge"        ); // charge						
  produces<vector<int> >	     ("mustrkcharge"	).setBranchAlias("mus_trk_charge"    ); // track charge						
  produces<vector<float> >	     ("musouterPhi"	).setBranchAlias("mus_outerPhi"      ); // phi angle of the outermost point in tracker		
  produces<vector<float> >	     ("musouterEta"	).setBranchAlias("mus_outerEta"      ); // eta angle of the outermost point in tracker		
  produces<vector<int> >             ("mustrkrefkey"    ).setBranchAlias("mus_trkrefkey"     ); // index of track from track ref stored in muon collection
  // muon quantities
  produces<vector<int> >             ("musnmatches"	).setBranchAlias("mus_nmatches"      ); // number of stations with matched segments                                                  
  produces<vector<float> >	     ("museem"		).setBranchAlias("mus_e_em"          ); // energy in crossed ECAL crystalls                                                          
  produces<vector<float> >	     ("musehad"		).setBranchAlias("mus_e_had"         ); // energy in crossed HCAL towers                                                             
  produces<vector<float> >	     ("museho"		).setBranchAlias("mus_e_ho"          ); // energy in crossed HO towers                                                               
  produces<vector<float> >	     ("museemS9"	).setBranchAlias("mus_e_emS9"        ); // energy in 3x3 ECAL crystall shape                                                         
  produces<vector<float> >	     ("musehadS9"	).setBranchAlias("mus_e_hadS9"       ); //energy in 3x3 HCAL towers                                                                 
  produces<vector<float> >	     ("musehoS9"	).setBranchAlias("mus_e_hoS9"        ); // energy in 3x3 HO towers                                                                   
  produces<vector<float> >           ("musiso"          ).setBranchAlias("mus_iso"           ); //mirrors the isolation in CMS1 (home grown trackIsolatio()
  produces<vector<float> >	     ("musiso03sumPt"	).setBranchAlias("mus_iso03_sumPt"   ); // sum of track Pt for cone of 0.3                                                           
  produces<vector<float> >	     ("musiso03emEt"	).setBranchAlias("mus_iso03_emEt"    ); // sum of ecal Et for cone of 0.3                                                            
  produces<vector<float> >	     ("musiso03hadEt"	).setBranchAlias("mus_iso03_hadEt"   ); // sum of hcal Et for cone of 0.3                                                            
  produces<vector<float> >	     ("musiso03hoEt"	).setBranchAlias("mus_iso03_hoEt"    ); // sum of ho Et for cone of 0.3                                                              
  produces<vector<int> >	     ("musiso03ntrk"	).setBranchAlias("mus_iso03_ntrk"    ); // number of tracks in the cone of 0.3                                                       
  produces<vector<float> >	     ("musiso05sumPt"	).setBranchAlias("mus_iso05_sumPt"   ); // sum of track Pt for cone of 0.5                                                           
  produces<vector<float> >	     ("musiso05emEt"    ).setBranchAlias("mus_iso05_emEt"    ); // sum of ecal Et for cone of 0.5                                                            
  produces<vector<float> >	     ("musiso05hadEt"	).setBranchAlias("mus_iso05_hadEt"   ); // sum of hcal Et for cone of 0.5                                                            
  produces<vector<float> >	     ("musiso05hoEt"	).setBranchAlias("mus_iso05_hoEt"    ); // sum of ho Et for cone of 0.5                                                              
  produces<vector<int> >	     ("musiso05ntrk"	).setBranchAlias("mus_iso05_ntrk"    ); // number of tracks in the cone of 0.5                                                       
        
  produces<vector<float> >	     ("musgfitchi2"     ).setBranchAlias("mus_gfit_chi2"     ); // chi2 of the global muon fit                                                               
  produces<vector<float> >	     ("musgfitndof"	).setBranchAlias("mus_gfit_ndof"     ); // number of degree of freedom of the global muon fit                                        
  produces<vector<int> >	     ("musgfitvalidHits").setBranchAlias("mus_gfit_validHits"); // number of valid hits of the global muon fit                
  produces<vector<Point> >            ("musgfitouterPos" ).setBranchAlias("mus_gfit_outerPos" ); //position of outermost hit
  // loose tracker muon identification based on muon/hadron penetration depth difference       
  produces<vector<int> >	     ("muspidTMLastStationLoose"    ).setBranchAlias("mus_pid_TMLastStationLoose"    ); 
  // tight tracker muon identification based on muon/hadron penetration depth difference       
  produces<vector<int> >	     ("muspidTMLastStationTight"    ).setBranchAlias("mus_pid_TMLastStationTight"    );
  // loose tracker muon likelihood identification based on muon matches and calo depositions   
  produces<vector<int> >	     ("muspidTM2DCompatibilityLoose").setBranchAlias("mus_pid_TM2DCompatibilityLoose");
  // tight tracker muon likelihood identification based on muon matches and calo depositions
  produces<vector<int> >	     ("muspidTM2DCompatibilityTight").setBranchAlias("mus_pid_TM2DCompatibilityTight"); 

  muonsInputTag  = iConfig.getParameter<edm::InputTag>("muonsInputTag" ); 
  tracksInputTag = iConfig.getParameter<edm::InputTag>("tracksInputTag");
  beamSpotInputTag = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
  
}

void MuonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  // make vectors to hold the information
  auto_ptr<vector<int> >	   vector_mus_type    	      (new vector<int>	           );        
  auto_ptr<vector<int> >	   vector_mus_goodmask       (new vector<int>             );        
  auto_ptr<vector<LorentzVector> > vector_mus_p4             (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > vector_mus_trk_p4	      (new vector<LorentzVector>   );
  auto_ptr<vector<float> >	   vector_mus_d0	      (new vector<float>           );      
  auto_ptr<vector<float> >	   vector_mus_z0	      (new vector<float>	   );      
  auto_ptr<vector<float> >	   vector_mus_d0corr	      (new vector<float>	   );      
  auto_ptr<vector<float> >	   vector_mus_z0corr	      (new vector<float>           );      
  auto_ptr<vector<float> >         vector_mus_vertexphi      (new vector<float>	   );      
  auto_ptr<vector<float> >	   vector_mus_chi2	      (new vector<float>	   );      
  auto_ptr<vector<float> >	   vector_mus_ndof	      (new vector<float>	   );      
  auto_ptr<vector<int> >	   vector_mus_validHits      (new vector<int>             );        
  auto_ptr<vector<int> >	   vector_mus_lostHits	      (new vector<int>             );        
  auto_ptr<vector<float> >	   vector_mus_d0Err	      (new vector<float>	   );      
  auto_ptr<vector<float> >	   vector_mus_z0Err	      (new vector<float>	   );      
  auto_ptr<vector<float> >	   vector_mus_ptErr	      (new vector<float>	   );      
  auto_ptr<vector<float> >	   vector_mus_etaErr	      (new vector<float>	   );      
  auto_ptr<vector<float> >	   vector_mus_phiErr	      (new vector<float>	   );      
  auto_ptr<vector<int> >	   vector_mus_charge	      (new vector<int>             );        
  auto_ptr<vector<int> >	   vector_mus_trk_charge     (new vector<int>             );        
  auto_ptr<vector<float> >	   vector_mus_outerPhi	      (new vector<float>	   );      
  auto_ptr<vector<float> >	   vector_mus_outerEta	      (new vector<float>	   );      
  auto_ptr<vector<int> >           vector_mus_trkrefkey      (new vector<int>             );
  auto_ptr<vector<int> >	   vector_mus_nmatches	      (new vector<int>             );
  auto_ptr<vector<float> >	   vector_mus_e_em	      (new vector<float>           );
  auto_ptr<vector<float> >	   vector_mus_e_had	      (new vector<float>           );
  auto_ptr<vector<float> >	   vector_mus_e_ho	      (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_e_emS9	      (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_e_hadS9	      (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_e_hoS9	      (new vector<float>	   );
  auto_ptr<vector<float> >         vector_mus_iso            (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_iso03_sumPt    (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_iso03_emEt     (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_iso03_hadEt    (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_iso03_hoEt     (new vector<float>	   );
  auto_ptr<vector<int> >	   vector_mus_iso03_ntrk     (new vector<int>  	   );
  auto_ptr<vector<float> >	   vector_mus_iso05_sumPt    (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_iso05_emEt     (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_iso05_hadEt    (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_iso05_hoEt     (new vector<float>	   );
  auto_ptr<vector<int> >	   vector_mus_iso05_ntrk     (new vector<int>  	   );
  auto_ptr<vector<float> >	   vector_mus_gfit_chi2      (new vector<float>	   );
  auto_ptr<vector<float> >	   vector_mus_gfit_ndof      (new vector<float>	   );
  auto_ptr<vector<int> >           vector_mus_gfit_validHits (new vector<int>  	   );
  auto_ptr<vector<Point> >         vector_mus_gfit_outerPos  (new vector<Point>           );
  auto_ptr<vector<int> >	   vector_mus_pid_TMLastStationLoose     (new vector<int> );
  auto_ptr<vector<int> >	   vector_mus_pid_TMLastStationTight     (new vector<int> );
  auto_ptr<vector<int> >	   vector_mus_pid_TM2DCompatibilityLoose (new vector<int> );
  auto_ptr<vector<int> >	   vector_mus_pid_TM2DCompatibilityTight (new vector<int> );

  // get muons
  Handle<edm::View<Muon> > muon_h;
  iEvent.getByLabel(muonsInputTag, muon_h);      // change this in the future
  edm::View<Muon>::const_iterator muons_end = muon_h->end();

  Handle<edm::View<Track> > tk_h;
  iEvent.getByLabel(tracksInputTag, tk_h);
  const edm::View<Track> *track_coll = tk_h.product();
     
  //get BeamSpot from BeamSpotMaker
  edm::InputTag beamSpot_tag(beamSpotInputTag.label(),"evtbs");
  edm::Handle<math::XYZPoint> beamSpotH;
  iEvent.getByLabel(beamSpot_tag, beamSpotH);
  const Point beamSpot = beamSpotH.isValid() ? *(beamSpotH.product()) : Point(0,0,0);

  for (edm::View<Muon>::const_iterator muon = muon_h->begin(); 
       muon != muons_end; ++muon) {
       
    float tempIso = trackRelIsolation(muon->momentum(), muon->vertex(), track_coll,
				      0.3,     //! dR < 0.3
				      0.01,    //! dR > 0.01
				      0.1,     //! |d0_tk| < 0.1 cm
				      999.9,   //! |el_2D - tk_2D| < 999
				      0.5,     //! |z0_el - z0_track| < 0.5
				      1.0,     //! min pt
				      7);      //! min nHits

       
    // fill vectors
    vector_mus_type         ->push_back(muon->type());
    int goodMask = 0;
    for (int iG = 0; iG < 16; ++iG){ //overkill here
      if (muon->isGood((Muon::SelectionType)iG) ) goodMask |= 
	(1 << iG);
    }
    
    const TrackRef siTrack     = muon->innerTrack();
    const TrackRef globalTrack = muon->globalTrack();

    vector_mus_goodmask      ->push_back(goodMask);
    vector_mus_p4            ->push_back(muon ->p4());
    vector_mus_trk_p4        ->push_back(siTrack.isNonnull() ? 
				 LorentzVector( siTrack.get()->px(), siTrack.get()->py(),
						siTrack.get()->pz(), siTrack.get()->p() )
				  : LorentzVector(0, 0, 0, 0));
    vector_mus_d0            ->push_back(siTrack.isNonnull() ? siTrack->d0()                            :  -999        );
    vector_mus_z0            ->push_back(siTrack.isNonnull() ? siTrack->dz()                            :  -999        );
    vector_mus_d0corr        ->push_back(siTrack.isNonnull() ? siTrack->dxy(beamSpot)                   :  -999        );
    vector_mus_z0corr        ->push_back(siTrack.isNonnull() ? siTrack->dxy(beamSpot)                   :  -999        );
    vector_mus_vertexphi     ->push_back(siTrack.isNonnull() ? atan2( siTrack->vy(), siTrack->vx() )    :  -999        );
    vector_mus_chi2          ->push_back(siTrack.isNonnull() ? siTrack->chi2()                          :  -999        );
    vector_mus_ndof          ->push_back(siTrack.isNonnull() ? siTrack->ndof()                          :  -999        );
    vector_mus_validHits     ->push_back(siTrack.isNonnull() ? siTrack->numberOfValidHits()             :  -999        );
    vector_mus_lostHits      ->push_back(siTrack.isNonnull() ? siTrack->numberOfLostHits()              :  -999        );
    vector_mus_d0Err         ->push_back(siTrack.isNonnull() ? siTrack->d0Error()                       :  -999        );
    vector_mus_z0Err         ->push_back(siTrack.isNonnull() ? siTrack->dzError()                       :  -999        );
    vector_mus_ptErr         ->push_back(siTrack.isNonnull() ? siTrack->ptError()                       :  -999        );
    vector_mus_etaErr        ->push_back(siTrack.isNonnull() ? siTrack->etaError()                      :  -999        );
    vector_mus_phiErr        ->push_back(siTrack.isNonnull() ? siTrack->phiError()                      :  -999        );
    vector_mus_charge        ->push_back(muon->charge()                                                                );
    vector_mus_trk_charge    ->push_back(siTrack.isNonnull() ? siTrack->charge()                        :  -999        );
    vector_mus_outerPhi      ->push_back(-999                                                                          );
    vector_mus_outerEta      ->push_back(-999                                                                          );
    vector_mus_trkrefkey     ->push_back(siTrack.isNonnull() ? (int)siTrack.index()                     :  -999        );
    vector_mus_nmatches      ->push_back(muon->isMatchesValid() ? muon->numberOfMatches()	         :  -999        );
    vector_mus_e_em          ->push_back(muon->isEnergyValid() ? muon->calEnergy().em   	         :  -999        );
    vector_mus_e_had         ->push_back(muon->isEnergyValid() ? muon->calEnergy().had		         :  -999        );
    vector_mus_e_ho          ->push_back(muon->isEnergyValid() ? muon->calEnergy().ho		         :  -999        );
    vector_mus_e_emS9        ->push_back(muon->isEnergyValid() ? muon->calEnergy().emS9		 :  -999        );
    vector_mus_e_hadS9       ->push_back(muon->isEnergyValid() ? muon->calEnergy().hadS9	         :  -999        );
    vector_mus_e_hoS9        ->push_back(muon->isEnergyValid() ? muon->calEnergy().hoS9                 :  -999        );
    vector_mus_iso           ->push_back(tempIso                                                                       );
    vector_mus_iso03_sumPt   ->push_back(muon->isIsolationValid() ? muon->isolationR03().sumPt           :  -999       );
    vector_mus_iso03_emEt    ->push_back(muon->isIsolationValid() ? muon->isolationR03().emEt	          :  -999       );
    vector_mus_iso03_hadEt   ->push_back(muon->isIsolationValid() ? muon->isolationR03().hadEt	          :  -999       );
    vector_mus_iso03_hoEt    ->push_back(muon->isIsolationValid() ? muon->isolationR03().hoEt	          :  -999       );
    vector_mus_iso03_ntrk    ->push_back(muon->isIsolationValid() ? muon->isolationR03().nTracks         :  -999       );
    vector_mus_iso05_sumPt   ->push_back(muon->isIsolationValid() ? muon->isolationR05().sumPt	          :  -999       );
    vector_mus_iso05_emEt    ->push_back(muon->isIsolationValid() ? muon->isolationR05().emEt	          :  -999       );
    vector_mus_iso05_hadEt   ->push_back(muon->isIsolationValid() ? muon->isolationR05().hadEt	          :  -999       );
    vector_mus_iso05_hoEt    ->push_back(muon->isIsolationValid() ? muon->isolationR05().hoEt	          :  -999       );
    vector_mus_iso05_ntrk    ->push_back(muon->isIsolationValid() ? muon->isolationR05().nTracks         :  -999       );
    vector_mus_gfit_chi2     ->push_back(globalTrack.isNonnull() ? muon->globalTrack()->chi2()	  :  -999       );
    vector_mus_gfit_ndof     ->push_back(globalTrack.isNonnull() ? muon->globalTrack()->ndof()	  :  -999       );
    vector_mus_gfit_validHits->push_back(globalTrack.isNonnull() ? muon->globalTrack()->numberOfValidHits() 	: -999	);
    vector_mus_gfit_outerPos ->push_back(globalTrack.isNonnull() ? muon->globalTrack()->outerPosition() :  Point(0,0,0));
    vector_mus_pid_TMLastStationLoose     ->push_back(muon->isMatchesValid() ? muon::isGoodMuon(*muon,Muon::TMLastStationLoose)     : -999	);
    vector_mus_pid_TMLastStationTight     ->push_back(muon->isMatchesValid() ? muon::isGoodMuon(*muon,Muon::TMLastStationTight)     : -999	);
    vector_mus_pid_TM2DCompatibilityLoose ->push_back(muon->isMatchesValid() ? muon::isGoodMuon(*muon,Muon::TM2DCompatibilityLoose)	: -999	);
    vector_mus_pid_TM2DCompatibilityTight ->push_back(muon->isMatchesValid() ? muon::isGoodMuon(*muon,Muon::TM2DCompatibilityTight)	: -999	);

  }
     
  // store vectors
  iEvent.put(vector_mus_type            , "mustype"              );
  iEvent.put(vector_mus_goodmask        , "musgoodmask"          );
  iEvent.put(vector_mus_p4              , "musp4"                );
  iEvent.put(vector_mus_trk_p4          , "mustrkp4"             );
  iEvent.put(vector_mus_d0              , "musd0"                );
  iEvent.put(vector_mus_z0              , "musz0"                );
  iEvent.put(vector_mus_d0corr          , "musd0corr"            );
  iEvent.put(vector_mus_z0corr          , "musz0corr"            );
  iEvent.put(vector_mus_vertexphi       , "musvertexphi"         );
  iEvent.put(vector_mus_chi2            , "muschi2"              );
  iEvent.put(vector_mus_ndof            , "musndof"              );
  iEvent.put(vector_mus_validHits       , "musvalidHits"         );
  iEvent.put(vector_mus_lostHits        , "muslostHits"          );
  iEvent.put(vector_mus_d0Err           , "musd0Err"             );
  iEvent.put(vector_mus_z0Err           , "musz0Err"             );
  iEvent.put(vector_mus_ptErr           , "musptErr"             );
  iEvent.put(vector_mus_etaErr          , "musetaErr"            );
  iEvent.put(vector_mus_phiErr          , "musphiErr"            );
  iEvent.put(vector_mus_charge          , "muscharge"            );
  iEvent.put(vector_mus_trk_charge      , "mustrkcharge"         );
  iEvent.put(vector_mus_outerPhi        , "musouterPhi"          );
  iEvent.put(vector_mus_outerEta        , "musouterEta"          );
  iEvent.put(vector_mus_trkrefkey       , "mustrkrefkey"         );
  iEvent.put(vector_mus_nmatches        , "musnmatches"          );
  iEvent.put(vector_mus_e_em            , "museem"               );
  iEvent.put(vector_mus_e_had		 , "musehad"              );
  iEvent.put(vector_mus_e_ho            , "museho"               );
  iEvent.put(vector_mus_e_emS9		 , "museemS9"             );
  iEvent.put(vector_mus_e_hadS9         , "musehadS9"            );
  iEvent.put(vector_mus_e_hoS9          , "musehoS9"             );
  iEvent.put(vector_mus_iso             , "musiso"               );
  iEvent.put(vector_mus_iso03_sumPt     , "musiso03sumPt"        );
  iEvent.put(vector_mus_iso03_emEt      , "musiso03emEt"         );
  iEvent.put(vector_mus_iso03_hadEt     , "musiso03hadEt"        );
  iEvent.put(vector_mus_iso03_hoEt      , "musiso03hoEt"         );
  iEvent.put(vector_mus_iso03_ntrk      , "musiso03ntrk"         );
  iEvent.put(vector_mus_iso05_sumPt     , "musiso05sumPt"        );
  iEvent.put(vector_mus_iso05_emEt      , "musiso05emEt"         );
  iEvent.put(vector_mus_iso05_hadEt     , "musiso05hadEt"        );
  iEvent.put(vector_mus_iso05_hoEt      , "musiso05hoEt"         );
  iEvent.put(vector_mus_iso05_ntrk      , "musiso05ntrk"         );
  iEvent.put(vector_mus_gfit_chi2       , "musgfitchi2"          );
  iEvent.put(vector_mus_gfit_ndof       , "musgfitndof"          );
  iEvent.put(vector_mus_gfit_validHits  , "musgfitvalidHits"     );
  iEvent.put(vector_mus_gfit_outerPos   , "musgfitouterPos"      );
  iEvent.put(vector_mus_pid_TMLastStationLoose	, "muspidTMLastStationLoose"     );
  iEvent.put(vector_mus_pid_TMLastStationTight	, "muspidTMLastStationTight"     );
  iEvent.put(vector_mus_pid_TM2DCompatibilityLoose	, "muspidTM2DCompatibilityLoose" );
  iEvent.put(vector_mus_pid_TM2DCompatibilityTight	, "muspidTM2DCompatibilityTight" );
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonMaker::endJob() {
}
//---------------------------------------------------------------------------
//Track Isolation
//---------------------------------------------------------------------------
double MuonMaker::trackRelIsolation(const math::XYZVector momentum,
				    const math::XYZPoint vertex,
				    const  edm::View<Track>* tracks,
				    double dRConeMax, double dRConeMin,
				    double tkVtxDMax,
				    double vtxDiffDMax, double vtxDiffZMax, double ptMin, unsigned int nHits)
{
  double isoResult = -10.;
  if ( tracks == 0 ) {
    cout << "Configuration Error: track collection is not set!" <<endl;
    return isoResult;
  }

  double sumPt = 0;

  edm::View<Track>::const_iterator iTk = tracks->begin();
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
DEFINE_FWK_MODULE(MuonMaker);
