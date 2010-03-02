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
// $Id: MuonMaker.cc,v 1.37 2010/03/02 19:36:08 fgolf Exp $
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

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Math/VectorUtil.h"

#include "DataFormats/Math/interface/Point3D.h"


typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace std;
using namespace reco;

MuonMaker::MuonMaker(const edm::ParameterSet& iConfig)
{
  // mu track quantities
  produces<vector<int> >	     ("mustype"		      ).setBranchAlias("mus_type"               );	// type
  produces<vector<int> >	     ("musgoodmask"	      ).setBranchAlias("mus_goodmask"           ); // good mask
  produces<vector<LorentzVector> >   ("musp4"		      ).setBranchAlias("mus_p4"                 ); // candidate p4->this can either be gfit p4, tracker p4 or STA p4 (only for STA muoons)						
  produces<vector<LorentzVector> >   ("mustrkp4"	      ).setBranchAlias("mus_trk_p4"             ); // track p4						
  produces<vector<LorentzVector> >   ("musgfitp4"             ).setBranchAlias("mus_gfit_p4"            ); // global fit p4, if global fit exists
  produces<vector<LorentzVector> >   ("musstap4"              ).setBranchAlias("mus_sta_p4"             ); // global fit p4, if global fit exists
  produces<vector<int>   >           ("mustrkidx"             ).setBranchAlias("mus_trkidx"             );	// track index matched to muon
  produces<vector<float> >	     ("musd0"		      ).setBranchAlias("mus_d0"                 ); // impact parameter at the point of closest approach	using the tracker fit
  produces<vector<float> >	     ("musz0"		      ).setBranchAlias("mus_z0"                 ); // z position of the point of closest approach. From the si track		
  produces<vector<float> >	     ("musd0corr"	      ).setBranchAlias("mus_d0corr"             ); // corrected impact parameter at the point of closest approach. From si track	
  produces<vector<float> >	     ("musz0corr"	      ).setBranchAlias("mus_z0corr"             ); // corrected z position of the point of closest approach. From si track		
  produces<vector<float> >	     ("musvertexphi"	      ).setBranchAlias("mus_vertexphi"          ); // phi angle of the point of closest approach. From si track		
  produces<vector<float> >	     ("muschi2"		      ).setBranchAlias("mus_chi2"               ); // chi2 of the silicon tracker fit			
  produces<vector<float> >	     ("musndof"		      ).setBranchAlias("mus_ndof"               ); // number of degrees of freedom of the si tracker fit		
  produces<vector<int> >	     ("musvalidHits"	      ).setBranchAlias("mus_validHits"          ); // number of used hits in the sitracker fit			
  produces<vector<int> >	     ("muslostHits"	      ).setBranchAlias("mus_lostHits"           ); // number of lost hits in the sitracker fit			
  produces<vector<int> >             ("musgfitvalidSTAHits"   ).setBranchAlias("mus_gfit_validSTAHits"  ); // number of hits in the stand alone fit that made it into the gfit
  produces<vector<int> >             ("musgfitvalidSiHits"    ).setBranchAlias("mus_gfit_validSiHits"   ); // number of hits in the Si fit that made it into the gfit
  produces<vector<float> >	     ("musd0Err"	      ).setBranchAlias("mus_d0Err"              ); // error on the impact parameter, si track fit			
  produces<vector<float> >	     ("musz0Err"	      ).setBranchAlias("mus_z0Err"              ); // error on z position of the point of closest approach, si track fit	
  produces<vector<float> >	     ("musptErr"	      ).setBranchAlias("mus_ptErr"              ); // si track Pt error					
  produces<vector<float> >	     ("musetaErr"	      ).setBranchAlias("mus_etaErr"             ); // si track eta error					
  produces<vector<float> >	     ("musphiErr"	      ).setBranchAlias("mus_phiErr"             ); // si track phi error					
  produces<vector<int> >	     ("muscharge"	      ).setBranchAlias("mus_charge"             ); // charge from muon object 						
  produces<vector<int> >	     ("mustrkcharge"	      ).setBranchAlias("mus_trk_charge"         ); // si track charge
  produces<vector<float> >           ("musqoverp"             ).setBranchAlias("mus_qoverp"             ); // si track qoverp
  produces<vector<float> >           ("musqoverpError"        ).setBranchAlias("mus_qoverpError"        ); // si track qoverp error
    // muon quantities
  produces<vector<int> >             ("musnmatches"	      ).setBranchAlias("mus_nmatches"           ); // number of stations with matched segments 
  produces<vector<float> >	     ("museem"		      ).setBranchAlias("mus_e_em"               ); // energy in crossed ECAL crystalls 
  produces<vector<float> >	     ("musehad"		      ).setBranchAlias("mus_e_had"              ); // energy in crossed HCAL towers 
  produces<vector<float> >	     ("museho"		      ).setBranchAlias("mus_e_ho"               ); // energy in crossed HO towers 
  produces<vector<float> >	     ("museemS9"	      ).setBranchAlias("mus_e_emS9"             ); // energy in 3x3 ECAL crystall shape 
  produces<vector<float> >	     ("musehadS9"	      ).setBranchAlias("mus_e_hadS9"            ); //energy in 3x3 HCAL towers 
  produces<vector<float> >	     ("musehoS9"	      ).setBranchAlias("mus_e_hoS9"             ); // energy in 3x3 HO towers 
  produces<vector<float> >           ("musisotrckvetoDep"     ).setBranchAlias("mus_iso_trckvetoDep"    );//sumPt in the veto cone, tracker
  produces<vector<float> >           ("musisoecalvetoDep"     ).setBranchAlias("mus_iso_ecalvetoDep"    );//sumEt in the veto cone, ecal
  produces<vector<float> >           ("musisohcalvetoDep"     ).setBranchAlias("mus_iso_hcalvetoDep"    );//sumPt in the veto cone, hcal
  produces<vector<float> >           ("musisohovetoDep"       ).setBranchAlias("mus_iso_hovetoDep"      );//sumPt in the veto cone, ho
  produces<vector<float> >	     ("musiso03sumPt"	      ).setBranchAlias("mus_iso03_sumPt"        ); // sum of track Pt for cone of 0.3 
  produces<vector<float> >	     ("musiso03emEt"	      ).setBranchAlias("mus_iso03_emEt"         ); // sum of ecal Et for cone of 0.3 
  produces<vector<float> >	     ("musiso03hadEt"	      ).setBranchAlias("mus_iso03_hadEt"        ); // sum of hcal Et for cone of 0.3 
  produces<vector<float> >	     ("musiso03hoEt"	      ).setBranchAlias("mus_iso03_hoEt"         ); // sum of ho Et for cone of 0.3 
  produces<vector<int> >	     ("musiso03ntrk"	      ).setBranchAlias("mus_iso03_ntrk"         ); // number of tracks in the cone of 0.3 
  produces<vector<float> >	     ("musiso05sumPt"	      ).setBranchAlias("mus_iso05_sumPt"        ); // sum of track Pt for cone of 0.5 
  produces<vector<float> >	     ("musiso05emEt"          ).setBranchAlias("mus_iso05_emEt"         ); // sum of ecal Et for cone of 0.5 
  produces<vector<float> >	     ("musiso05hadEt"	      ).setBranchAlias("mus_iso05_hadEt"        ); // sum of hcal Et for cone of 0.5 
  produces<vector<float> >	     ("musiso05hoEt"	      ).setBranchAlias("mus_iso05_hoEt"         ); // sum of ho Et for cone of 0.5 
  produces<vector<int> >	     ("musiso05ntrk"	      ).setBranchAlias("mus_iso05_ntrk"         ); // number of tracks in the cone of 0.5 

  //new
  produces<vector<float> >           ("musgfitd0"             ).setBranchAlias("mus_gfit_d0"            ); // d0 from global fit, if it exists
  produces<vector<float> >           ("musgfitz0"             ).setBranchAlias("mus_gfit_z0"            ); // z0 from global fit, if it exists
  produces<vector<float> >           ("musgfitd0Err"          ).setBranchAlias("mus_gfit_d0Err"         ); // d0Err from global fit, if it exists
  produces<vector<float> >           ("musgfitz0Err"          ).setBranchAlias("mus_gfit_z0Err"         ); // z0Err from global fit, if it exists
  produces<vector<float> >           ("musgfitd0corr"         ).setBranchAlias("mus_gfit_d0corr"        ); // Beamspot corrected d0 from global fit, if it exists
  produces<vector<float> >           ("musgfitz0corr"         ).setBranchAlias("mus_gfit_z0corr"        ); // Beamspot corrected z0 from global fit, if it exists
  produces<vector<float> >           ("musgfitqoverp"         ).setBranchAlias("mus_gfit_qoverp"        ); // global track qoverp
  produces<vector<float> >           ("musgfitqoverpError"    ).setBranchAlias("mus_gfit_qoverpError"   ); // global track qoverp error  
  
  produces<vector<float> >	     ("musgfitchi2"           ).setBranchAlias("mus_gfit_chi2"          ); // chi2 of the global muon fit 
  produces<vector<float> >	     ("musgfitndof"	      ).setBranchAlias("mus_gfit_ndof"          ); // number of degree of freedom of the global muon fit 
  produces<vector<int> >	     ("musgfitvalidHits"      ).setBranchAlias("mus_gfit_validHits"     ); // number of valid hits of the global muon fit 
  //STA fit crap
  produces<vector<float> >           ("musstad0"             ).setBranchAlias("mus_sta_d0"            ); // d0 from STA fit, if it exists
  produces<vector<float> >           ("musstaz0"             ).setBranchAlias("mus_sta_z0"            ); // z0 from STA fit, if it exists
  produces<vector<float> >           ("musstad0Err"          ).setBranchAlias("mus_sta_d0Err"         ); // d0Err from STA fit, if it exists
  produces<vector<float> >           ("musstaz0Err"          ).setBranchAlias("mus_sta_z0Err"         ); // z0Err from STA fit, if it exists
  produces<vector<float> >           ("musstad0corr"         ).setBranchAlias("mus_sta_d0corr"        ); // Beamspot corrected d0 from STA fit, if it exists
  produces<vector<float> >           ("musstaz0corr"         ).setBranchAlias("mus_sta_z0corr"        ); // Beamspot corrected z0 from STA fit, if it exists
  produces<vector<float> >           ("musstaqoverp"         ).setBranchAlias("mus_sta_qoverp"        ); // STA track qoverp
  produces<vector<float> >           ("musstaqoverpError"    ).setBranchAlias("mus_sta_qoverpError"   ); // STA track qoverp error  
  
  produces<vector<float> >	     ("musstachi2"           ).setBranchAlias("mus_sta_chi2"          ); // chi2 of the STA muon fit 
  produces<vector<float> >	     ("musstandof"	     ).setBranchAlias("mus_sta_ndof"          ); // number of degree of freedom of the STA muon fit 
  produces<vector<int> >	     ("musstavalidHits"      ).setBranchAlias("mus_sta_validHits"     ); // number of valid hits of the STA muon fit 
  

  
  //Muon timing info -> http://cmslxr.fnal.gov/lxr/source/DataFormats/MuonReco/interface/MuonTime.h
  produces<vector<int> >             ("mustimeNumStationsUsed").setBranchAlias("mus_timeNumStationsUsed"); //number of muon stations used for timing info
  // time of arrival at the IP for the Beta=1 hypothesis -> particle moving from inside out
  produces<vector<float> >           ("mustimeAtIpInOut"      ).setBranchAlias("mus_timeAtIpInOut"      ); 
  produces<vector<float> >           ("mustimeAtIpInOutErr"   ).setBranchAlias("mus_timeAtIpInOutErr"   );
  //particle moving from outside in
  produces<vector<float> >           ("mustimeAtIpOutIn"      ).setBranchAlias("mus_timeAtIpOutIn"      );
  produces<vector<float> >           ("mustimeAtIpOutInErr"   ).setBranchAlias("mus_timeAtIpOutInErr"   );
  //direction estimate based on time dispersion. Enum defn given in above header file. in 312: 
  //Direction { OutsideIn = -1, Undefined = 0, InsideOut = 1 };
  produces<vector<int> >             ("mustimeDirection"      ).setBranchAlias("mus_timeDirection"      );

  
  
  // loose tracker muon identification based on muon/hadron penetration depth difference       
  produces<vector<int> >	     ("muspidTMLastStationLoose"    ).setBranchAlias("mus_pid_TMLastStationLoose"    ); 
  // tight tracker muon identification based on muon/hadron penetration depth difference       
  produces<vector<int> >	     ("muspidTMLastStationTight"    ).setBranchAlias("mus_pid_TMLastStationTight"    );
  // loose tracker muon likelihood identification based on muon matches and calo depositions   
  produces<vector<int> >	     ("muspidTM2DCompatibilityLoose").setBranchAlias("mus_pid_TM2DCompatibilityLoose");
  // tight tracker muon likelihood identification based on muon matches and calo depositions
  produces<vector<int> >	     ("muspidTM2DCompatibilityTight").setBranchAlias("mus_pid_TM2DCompatibilityTight"); 
  //calo compatibility variable
  produces<vector<float> >           ("muscaloCompatibility").setBranchAlias("mus_caloCompatibility");
  //p4 because we're not able to (yet) read XYZPointDs in bare root for some reason 
  //the 4th co-ordinate is 0
  produces<vector<LorentzVector> >   ("musvertexp4"         ).setBranchAlias("mus_vertex_p4" ); // from the silicon fit
  produces<vector<LorentzVector> >   ("musgfitvertexp4"     ).setBranchAlias("mus_gfit_vertex_p4");
  produces<vector<LorentzVector> >   ("musgfitouterPosp4"   ).setBranchAlias("mus_gfit_outerPos_p4");
  produces<vector<LorentzVector> >   ("musstavertexp4"      ).setBranchAlias("mus_sta_vertex_p4");
  produces<vector<LorentzVector> >   ("musfitdefaultp4" ).setBranchAlias("mus_fitdefault_p4" );
  produces<vector<LorentzVector> >   ("musfitfirsthitp4").setBranchAlias("mus_fitfirsthit_p4");
  produces<vector<LorentzVector> >   ("musfitpickyp4"   ).setBranchAlias("mus_fitpicky_p4"	 );
  produces<vector<LorentzVector> >   ("musfittevp4"     ).setBranchAlias("mus_fittev_p4"     );
  
  muonsInputTag  		= iConfig.getParameter<edm::InputTag>("muonsInputTag" ); 
  beamSpotInputTag 		= iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
  tevMuonsName    		= iConfig.getParameter<string>("tevMuonsName" ); 
 
}

void MuonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  // make vectors to hold the information
  auto_ptr<vector<int> >	   vector_mus_type    	          (new vector<int>	       );        
  auto_ptr<vector<int> >	   vector_mus_goodmask            (new vector<int>             );        
  auto_ptr<vector<LorentzVector> > vector_mus_p4                  (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > vector_mus_trk_p4	          (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > vector_mus_gfit_p4             (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > vector_mus_sta_p4              (new vector<LorentzVector>   );
  auto_ptr<vector<int>   >         vector_mus_trkidx              (new vector<int>             );
  auto_ptr<vector<float> >	   vector_mus_d0	          (new vector<float>           );      
  auto_ptr<vector<float> >	   vector_mus_z0	          (new vector<float>	       );      
  auto_ptr<vector<float> >	   vector_mus_d0corr	          (new vector<float>	       );      
  auto_ptr<vector<float> >	   vector_mus_z0corr	          (new vector<float>           );      
  auto_ptr<vector<float> >         vector_mus_vertexphi           (new vector<float>	       );      
  auto_ptr<vector<float> >	   vector_mus_chi2	          (new vector<float>	       );      
  auto_ptr<vector<float> >	   vector_mus_ndof	          (new vector<float>	       );      
  auto_ptr<vector<int> >	   vector_mus_validHits           (new vector<int>             );        
  auto_ptr<vector<int> >	   vector_mus_lostHits	          (new vector<int>             );        
  auto_ptr<vector<int> >           vector_mus_gfit_validSTAHits   (new vector<int>             );
  auto_ptr<vector<int> >           vector_mus_gfit_validSiHits    (new vector<int>             );
  auto_ptr<vector<float> >	   vector_mus_d0Err	          (new vector<float>           );      
  auto_ptr<vector<float> >	   vector_mus_z0Err	          (new vector<float>	       );      
  auto_ptr<vector<float> >	   vector_mus_ptErr	          (new vector<float>	       );      
  auto_ptr<vector<float> >	   vector_mus_etaErr	          (new vector<float>	       );      
  auto_ptr<vector<float> >	   vector_mus_phiErr	          (new vector<float>	       );      
  auto_ptr<vector<int> >	   vector_mus_charge	          (new vector<int>             );        
  auto_ptr<vector<int> >	   vector_mus_trk_charge          (new vector<int>             );   
  auto_ptr<vector<float> >         vector_mus_qoverp              (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_qoverpError         (new vector<float>           );
  auto_ptr<vector<int> >	   vector_mus_nmatches	          (new vector<int>             );
  auto_ptr<vector<float> >	   vector_mus_e_em	          (new vector<float>           );
  auto_ptr<vector<float> >	   vector_mus_e_had   	          (new vector<float>           );
  auto_ptr<vector<float> >	   vector_mus_e_ho	          (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_e_emS9	          (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_e_hadS9	          (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_e_hoS9	          (new vector<float>	       );
  auto_ptr<vector<float> >         vector_mus_iso_trckvetoDep     (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_iso_ecalvetoDep     (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_iso_hcalvetoDep     (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_iso_hovetoDep       (new vector<float>           );
  auto_ptr<vector<float> >	   vector_mus_iso03_sumPt         (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_iso03_emEt          (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_iso03_hadEt         (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_iso03_hoEt          (new vector<float>	       );
  auto_ptr<vector<int> >	   vector_mus_iso03_ntrk          (new vector<int>  	       );
  auto_ptr<vector<float> >	   vector_mus_iso05_sumPt         (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_iso05_emEt          (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_iso05_hadEt         (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_iso05_hoEt          (new vector<float>	       );
  auto_ptr<vector<int> >	   vector_mus_iso05_ntrk          (new vector<int>  	       );
  //gfit
  auto_ptr<vector<float> >         vector_mus_gfit_d0             (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_gfit_z0             (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_gfit_d0Err          (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_gfit_z0Err          (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_gfit_d0corr         (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_gfit_z0corr         (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_gfit_qoverp         (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_gfit_qoverpError    (new vector<float>           );
  auto_ptr<vector<float> >	   vector_mus_gfit_chi2           (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_gfit_ndof           (new vector<float>	       );
  auto_ptr<vector<int> >           vector_mus_gfit_validHits      (new vector<int>  	       );
  //sta
  auto_ptr<vector<float> >         vector_mus_sta_d0             (new vector<float>            );
  auto_ptr<vector<float> >         vector_mus_sta_z0             (new vector<float>            );
  auto_ptr<vector<float> >         vector_mus_sta_d0Err          (new vector<float>            );
  auto_ptr<vector<float> >         vector_mus_sta_z0Err          (new vector<float>            );
  auto_ptr<vector<float> >         vector_mus_sta_d0corr         (new vector<float>            );
  auto_ptr<vector<float> >         vector_mus_sta_z0corr         (new vector<float>            );
  auto_ptr<vector<float> >         vector_mus_sta_qoverp         (new vector<float>            );
  auto_ptr<vector<float> >         vector_mus_sta_qoverpError    (new vector<float>            );
  auto_ptr<vector<float> >	   vector_mus_sta_chi2           (new vector<float>	       );
  auto_ptr<vector<float> >	   vector_mus_sta_ndof           (new vector<float>	       );
  auto_ptr<vector<int> >           vector_mus_sta_validHits      (new vector<int>  	       );

  
  auto_ptr<vector<int> >           vector_mus_timeNumStationsUsed (new vector<int>             );
  auto_ptr<vector<float> >         vector_mus_timeAtIpInOut       (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_timeAtIpInOutErr    (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_timeAtIpOutIn       (new vector<float>           );
  auto_ptr<vector<float> >         vector_mus_timeAtIpOutInErr    (new vector<float>           );
  auto_ptr<vector<int> >           vector_mus_timeDirection       (new vector<int>             );
  auto_ptr<vector<int> >	   vector_mus_pid_TMLastStationLoose       (new vector<int>    );
  auto_ptr<vector<int> >	   vector_mus_pid_TMLastStationTight       (new vector<int>    );
  auto_ptr<vector<int> >	   vector_mus_pid_TM2DCompatibilityLoose   (new vector<int>    );
  auto_ptr<vector<int> >	   vector_mus_pid_TM2DCompatibilityTight   (new vector<int>    );
  auto_ptr<vector<float> >         vector_mus_caloCompatibility            (new vector<float>  );

  auto_ptr<vector<LorentzVector> > vector_mus_vertex_p4                    (new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_gfit_vertex_p4               (new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_gfit_outerPos_p4             (new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_sta_vertex_p4                (new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_fitdefault_p4                (new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_fitfirsthit_p4               (new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_fitpicky_p4                  (new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_fittev_p4                    (new vector<LorentzVector> );
    
  // get muons
  Handle<edm::View<Muon> > muon_h;
  iEvent.getByLabel(muonsInputTag, muon_h);      // change this in the future
  edm::View<Muon>::const_iterator muons_end = muon_h->end();

  //get BeamSpot from BeamSpotMaker
  edm::InputTag beamSpot_tag(beamSpotInputTag.label(),"evtbsp4");
  edm::Handle<LorentzVector> beamSpotH;
  iEvent.getByLabel(beamSpot_tag, beamSpotH);
  const Point beamSpot = beamSpotH.isValid() ?
                         Point(beamSpotH->x(), beamSpotH->y(), beamSpotH->z()) : Point(0,0,0);

  //maps for alternative muon fits
  Handle<TrackToTrackMap> trackMap;
  Handle<TrackToTrackMap> trackMapDefault;
  Handle<TrackToTrackMap> trackMapFirstHit;
  Handle<TrackToTrackMap> trackMapPicky;
  iEvent.getByLabel(tevMuonsName, "default", trackMapDefault);
  iEvent.getByLabel(tevMuonsName, "firstHit", trackMapFirstHit);
  iEvent.getByLabel(tevMuonsName, "picky", trackMapPicky);
  
  for (edm::View<Muon>::const_iterator muon = muon_h->begin(); 
       muon != muons_end; ++muon) {
    
    const TrackRef siTrack     = muon->innerTrack();
    const TrackRef globalTrack = muon->globalTrack();
    const TrackRef staTrack    = muon->outerTrack();

    // fill vectors
    vector_mus_type         ->push_back(muon->type());
    int goodMask = 0;
    
    for (int iG = 0; iG < 16; ++iG){ //overkill here
      if (isGoodMuon(*muon,(muon::SelectionType)iG) ) goodMask |=   (1 << iG);
    }
     
    vector_mus_goodmask           ->push_back(goodMask);
    vector_mus_p4                 ->push_back( LorentzVector( muon   ->p4() ) );
    vector_mus_trk_p4             ->push_back(siTrack.isNonnull() ? 
				 LorentzVector( siTrack.get()->px(), siTrack.get()->py(),
						siTrack.get()->pz(), siTrack.get()->p() )
				  : LorentzVector(0, 0, 0, 0));
    vector_mus_gfit_p4            ->push_back( globalTrack.isNonnull() ? 
					       LorentzVector(globalTrack->px(), globalTrack->py(),
							     globalTrack->pz(), globalTrack->p()) 
					       : LorentzVector(0, 0, 0, 0) );
    vector_mus_sta_p4             ->push_back( staTrack.isNonnull() ? 
					       LorentzVector(staTrack->px(), staTrack->py(),
							     staTrack->pz(), staTrack->p()) 
					       : LorentzVector(0, 0, 0, 0) );
					       
    vector_mus_trkidx             ->push_back(siTrack.isNonnull() ? static_cast<int>(siTrack.key())          :  -9999        );
    vector_mus_d0                 ->push_back(siTrack.isNonnull() ? siTrack->d0()                            :  -9999.       );
    vector_mus_z0                 ->push_back(siTrack.isNonnull() ? siTrack->dz()                            :  -9999.       );
    vector_mus_d0corr             ->push_back(siTrack.isNonnull() ? -1*(siTrack->dxy(beamSpot))              :  -9999.       );
    vector_mus_z0corr             ->push_back(siTrack.isNonnull() ? siTrack->dz(beamSpot)                    :  -9999.       );
    vector_mus_vertexphi          ->push_back(siTrack.isNonnull() ? atan2( siTrack->vy(), siTrack->vx() )    :  -9999.       );
    vector_mus_chi2               ->push_back(siTrack.isNonnull() ? siTrack->chi2()                          :  -9999.       );
    vector_mus_ndof               ->push_back(siTrack.isNonnull() ? siTrack->ndof()                          :  -9999.       );
    vector_mus_validHits          ->push_back(siTrack.isNonnull() ? siTrack->numberOfValidHits()             :  -9999        );
    vector_mus_lostHits           ->push_back(siTrack.isNonnull() ? siTrack->numberOfLostHits()              :  -9999        );
    vector_mus_gfit_validSTAHits  ->push_back(globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfValidMuonHits()   : -9999);
    vector_mus_gfit_validSiHits  ->push_back(globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfValidTrackerHits(): -9999);
    vector_mus_d0Err              ->push_back(siTrack.isNonnull() ? siTrack->d0Error()                       :  -9999.        );
    vector_mus_z0Err              ->push_back(siTrack.isNonnull() ? siTrack->dzError()                       :  -9999.        );
    vector_mus_ptErr              ->push_back(siTrack.isNonnull() ? siTrack->ptError()                       :  -9999.        );
    vector_mus_etaErr             ->push_back(siTrack.isNonnull() ? siTrack->etaError()                      :  -9999.        );
    vector_mus_phiErr             ->push_back(siTrack.isNonnull() ? siTrack->phiError()                      :  -9999.        );
    vector_mus_charge             ->push_back(muon->charge()                                                                  );
    vector_mus_trk_charge         ->push_back(siTrack.isNonnull() ? siTrack->charge()                        :  -9999         );
    vector_mus_qoverp             ->push_back(siTrack.isNonnull() ? siTrack->qoverp()                        :  -9999.        );
    vector_mus_qoverpError        ->push_back(siTrack.isNonnull() ? siTrack->qoverpError()                   :  -9999.        );
    vector_mus_nmatches           ->push_back(muon->isMatchesValid() ? muon->numberOfMatches()               :  -9999         );
    vector_mus_e_em               ->push_back(muon->isEnergyValid() ? muon->calEnergy().em                   :  -9999.        );
    vector_mus_e_had              ->push_back(muon->isEnergyValid() ? muon->calEnergy().had		     :  -9999.        );
    vector_mus_e_ho               ->push_back(muon->isEnergyValid() ? muon->calEnergy().ho		     :  -9999.        );
    vector_mus_e_emS9             ->push_back(muon->isEnergyValid() ? muon->calEnergy().emS9		     :  -9999.        );
    vector_mus_e_hadS9            ->push_back(muon->isEnergyValid() ? muon->calEnergy().hadS9	             :  -9999.        );
    vector_mus_e_hoS9             ->push_back(muon->isEnergyValid() ? muon->calEnergy().hoS9                 :  -9999.        );
    vector_mus_iso_trckvetoDep    ->push_back(muon->isEnergyValid() ? muon->isolationR03().trackerVetoPt     :  -9999.        );
    vector_mus_iso_ecalvetoDep    ->push_back(muon->isEnergyValid() ? muon->isolationR03().emVetoEt          :  -9999.        );      
    vector_mus_iso_hcalvetoDep    ->push_back(muon->isEnergyValid() ? muon->isolationR03().hadVetoEt         :  -9999.        );      
    vector_mus_iso_hovetoDep      ->push_back(muon->isEnergyValid() ? muon->isolationR03().hoVetoEt          :  -9999.        );      
    
    vector_mus_iso03_sumPt        ->push_back(muon->isIsolationValid() ? muon->isolationR03().sumPt          :  -9999.        );
    vector_mus_iso03_emEt         ->push_back(muon->isIsolationValid() ? muon->isolationR03().emEt	     :  -9999.        );
    vector_mus_iso03_hadEt        ->push_back(muon->isIsolationValid() ? muon->isolationR03().hadEt	     :  -9999.        );
    vector_mus_iso03_hoEt         ->push_back(muon->isIsolationValid() ? muon->isolationR03().hoEt	     :  -9999.        );
    vector_mus_iso03_ntrk         ->push_back(muon->isIsolationValid() ? muon->isolationR03().nTracks        :  -9999         );
    vector_mus_iso05_sumPt        ->push_back(muon->isIsolationValid() ? muon->isolationR05().sumPt          :  -9999.        );
    vector_mus_iso05_emEt         ->push_back(muon->isIsolationValid() ? muon->isolationR05().emEt	     :  -9999.        );
    vector_mus_iso05_hadEt        ->push_back(muon->isIsolationValid() ? muon->isolationR05().hadEt	     :  -9999.        );
    vector_mus_iso05_hoEt         ->push_back(muon->isIsolationValid() ? muon->isolationR05().hoEt	     :  -9999.        );
    vector_mus_iso05_ntrk         ->push_back(muon->isIsolationValid() ? muon->isolationR05().nTracks        :  -9999         );

    vector_mus_gfit_d0            ->push_back(globalTrack.isNonnull()  ? globalTrack->d0()                   :  -9999.        );
    vector_mus_gfit_z0            ->push_back(globalTrack.isNonnull()  ? globalTrack->dz()                   :  -9999.        );
    vector_mus_gfit_d0Err         ->push_back(globalTrack.isNonnull()  ? globalTrack->d0Error()              :  -9999.        );
    vector_mus_gfit_z0Err         ->push_back(globalTrack.isNonnull()  ? globalTrack->dzError()              :  -9999.        );
    vector_mus_gfit_d0corr        ->push_back(globalTrack.isNonnull()  ? -1*(globalTrack->dxy(beamSpot))     :  -9999.        );
    vector_mus_gfit_z0corr        ->push_back(globalTrack.isNonnull()  ? globalTrack->dz(beamSpot)           :  -9999.        );
    vector_mus_gfit_qoverp        ->push_back(globalTrack.isNonnull()  ? globalTrack->qoverp()               :  -9999.        );
    vector_mus_gfit_qoverpError   ->push_back(globalTrack.isNonnull()  ? globalTrack->qoverpError()          :  -9999.        );
    vector_mus_gfit_chi2          ->push_back(globalTrack.isNonnull()  ? globalTrack->chi2()	             :  -9999.        );
    vector_mus_gfit_ndof          ->push_back(globalTrack.isNonnull()  ? globalTrack->ndof()	             :  -9999         );
    vector_mus_gfit_validHits     ->push_back(globalTrack.isNonnull()  ? globalTrack->numberOfValidHits()    :  -9999         );

    //STA crap
    vector_mus_sta_d0            ->push_back(staTrack.isNonnull()  ? staTrack->d0()                   :  -9999.        );
    vector_mus_sta_z0            ->push_back(staTrack.isNonnull()  ? staTrack->dz()                   :  -9999.        );
    vector_mus_sta_d0Err         ->push_back(staTrack.isNonnull()  ? staTrack->d0Error()              :  -9999.        );
    vector_mus_sta_z0Err         ->push_back(staTrack.isNonnull()  ? staTrack->dzError()              :  -9999.        );
    vector_mus_sta_d0corr        ->push_back(staTrack.isNonnull()  ? -1*(staTrack->dxy(beamSpot))     :  -9999.        );
    vector_mus_sta_z0corr        ->push_back(staTrack.isNonnull()  ? staTrack->dz(beamSpot)           :  -9999.        );
    vector_mus_sta_qoverp        ->push_back(staTrack.isNonnull()  ? staTrack->qoverp()               :  -9999.        );
    vector_mus_sta_qoverpError   ->push_back(staTrack.isNonnull()  ? staTrack->qoverpError()          :  -9999.        );
    vector_mus_sta_chi2          ->push_back(staTrack.isNonnull()  ? staTrack->chi2()	            :  -9999.        );
    vector_mus_sta_ndof          ->push_back(staTrack.isNonnull()  ? staTrack->ndof()	            :  -9999         );
    vector_mus_sta_validHits     ->push_back(staTrack.isNonnull()  ? staTrack->numberOfValidHits()    :  -9999         );

    bool timeIsValid = muon->isTimeValid();
    vector_mus_timeNumStationsUsed->push_back(timeIsValid              ?  muon->time().nDof                   :  -9999         );
    vector_mus_timeAtIpInOut      ->push_back(timeIsValid              ?  muon->time().timeAtIpInOut          :  -9999.        );
    vector_mus_timeAtIpInOutErr   ->push_back(timeIsValid              ?  muon->time().timeAtIpInOutErr       :  -9999.        );
    vector_mus_timeAtIpOutIn      ->push_back(timeIsValid              ?  muon->time().timeAtIpOutIn          :  -9999.        );
    vector_mus_timeAtIpOutInErr   ->push_back(timeIsValid              ?  muon->time().timeAtIpOutInErr       :  -9999.        );
    vector_mus_timeDirection      ->push_back(timeIsValid              ?  muon->time().direction()            :  -9999         );
    vector_mus_pid_TMLastStationLoose     ->push_back(muon->isMatchesValid() ? muon::isGoodMuon(*muon,muon::TMLastStationLoose)    : -9999	);
    vector_mus_pid_TMLastStationTight     ->push_back(muon->isMatchesValid() ? muon::isGoodMuon(*muon,muon::TMLastStationTight)    : -9999	);
    vector_mus_pid_TM2DCompatibilityLoose ->push_back(muon->isMatchesValid() ? muon::isGoodMuon(*muon,muon::TM2DCompatibilityLoose): -9999	);
    vector_mus_pid_TM2DCompatibilityTight ->push_back(muon->isMatchesValid() ? muon::isGoodMuon(*muon,muon::TM2DCompatibilityTight): -9999	);
    vector_mus_caloCompatibility          ->push_back(muon->caloCompatibility() );

    vector_mus_vertex_p4                  ->push_back(siTrack.isNonnull() ? 
						      LorentzVector(siTrack->vx(),
								    siTrack->vy(),
								    siTrack->vz(), 0.) 
						      : LorentzVector(-9999.,-9999.,-9999.,-9999.) );
    vector_mus_gfit_vertex_p4             ->push_back(globalTrack.isNonnull() ? 
						      LorentzVector(globalTrack->vx(),
								    globalTrack->vy(),
								    globalTrack->vz(), 0.) 
						      : LorentzVector(-9999.,-9999.,-9999.,-9999.) );
    vector_mus_gfit_outerPos_p4          ->push_back(globalTrack.isNonnull() ?
						     LorentzVector(globalTrack->outerPosition().x(),
								   globalTrack->outerPosition().y(),
								   globalTrack->outerPosition().z(),0. )
						     : LorentzVector(-9999.,-9999.,-9999.,-9999.) );
    
	// if muon is not global
    if( !muon->isGlobalMuon() ) {
      vector_mus_fitdefault_p4 ->push_back( LorentzVector( 0, 0, 0, 0 ) );
      vector_mus_fitfirsthit_p4->push_back( LorentzVector( 0, 0, 0, 0 ) );
      vector_mus_fitpicky_p4   ->push_back( LorentzVector( 0, 0, 0, 0 ) );
      vector_mus_fittev_p4     ->push_back( LorentzVector( 0, 0, 0, 0 ) );
    }

	// if muon is global
    else {

      reco::TrackToTrackMap::const_iterator fittmp;

      if( !muon->combinedMuon().isAvailable() )
	std::cout << "WTF" << std::endl;

      fittmp = (*trackMapDefault).find(muon->combinedMuon());
      if( fittmp != trackMapDefault->end()  )
	vector_mus_fitdefault_p4->push_back( LorentzVector( (*fittmp).val->px(), (*fittmp).val->py(), (*fittmp).val->pz(), (*fittmp).val->p() ) );
      else
	vector_mus_fitdefault_p4 ->push_back( LorentzVector( 0, 0, 0, 0 ) );

      fittmp = (*trackMapFirstHit).find(muon->combinedMuon());
      if( fittmp != trackMapFirstHit->end()  )
	vector_mus_fitfirsthit_p4->push_back( LorentzVector( (*fittmp).val->px(), (*fittmp).val->py(), (*fittmp).val->pz(), (*fittmp).val->p() ) );
      else
	vector_mus_fitfirsthit_p4->push_back( LorentzVector( 0, 0, 0, 0 ) );

      fittmp = (*trackMapPicky).find(muon->combinedMuon());
      if( fittmp != trackMapPicky->end()  )
	vector_mus_fitpicky_p4->push_back( LorentzVector( (*fittmp).val->px(), (*fittmp).val->py(), (*fittmp).val->pz(), (*fittmp).val->p() ) );
      else
	vector_mus_fitpicky_p4->push_back( LorentzVector( 0, 0, 0, 0 ) );
      
      TrackRef fittmpref;

      fittmpref = muon::tevOptimized(*muon, *trackMapDefault, *trackMapFirstHit, *trackMapPicky);

      if( fittmpref.isAvailable() )
	vector_mus_fittev_p4->push_back( LorentzVector( fittmpref->px(), fittmpref->py(), fittmpref->pz(), fittmpref->p() ) );
      else
	vector_mus_fittev_p4     ->push_back( LorentzVector( 0, 0, 0, 0 ) );
    }
  }
     
  // store vectors
  iEvent.put(vector_mus_type               , "mustype"               );
  iEvent.put(vector_mus_goodmask           , "musgoodmask"           );
  iEvent.put(vector_mus_p4                 , "musp4"                 );
  iEvent.put(vector_mus_trk_p4             , "mustrkp4"              );
  iEvent.put(vector_mus_gfit_p4            , "musgfitp4"             );
  iEvent.put(vector_mus_sta_p4             , "musstap4"              );
  iEvent.put(vector_mus_trkidx             , "mustrkidx"             );
  iEvent.put(vector_mus_d0                 , "musd0"                 );
  iEvent.put(vector_mus_z0                 , "musz0"                 );
  iEvent.put(vector_mus_d0corr             , "musd0corr"             );
  iEvent.put(vector_mus_z0corr             , "musz0corr"             );
  iEvent.put(vector_mus_vertexphi          , "musvertexphi"          );
  iEvent.put(vector_mus_chi2               , "muschi2"               );
  iEvent.put(vector_mus_ndof               , "musndof"               );
  iEvent.put(vector_mus_validHits          , "musvalidHits"          );
  iEvent.put(vector_mus_lostHits           , "muslostHits"           );
  iEvent.put(vector_mus_gfit_validSTAHits  , "musgfitvalidSTAHits"   );
  iEvent.put(vector_mus_gfit_validSiHits   , "musgfitvalidSiHits"    );
  iEvent.put(vector_mus_d0Err              , "musd0Err"              );
  iEvent.put(vector_mus_z0Err              , "musz0Err"              );
  iEvent.put(vector_mus_ptErr              , "musptErr"              );
  iEvent.put(vector_mus_etaErr             , "musetaErr"             );
  iEvent.put(vector_mus_phiErr             , "musphiErr"             );
  iEvent.put(vector_mus_charge             , "muscharge"             );
  iEvent.put(vector_mus_trk_charge         , "mustrkcharge"          );
  iEvent.put(vector_mus_qoverp             , "musqoverp"             );
  iEvent.put(vector_mus_qoverpError        , "musqoverpError"        );
  iEvent.put(vector_mus_nmatches           , "musnmatches"           );
  iEvent.put(vector_mus_e_em               , "museem"                );
  iEvent.put(vector_mus_e_had              , "musehad"               );
  iEvent.put(vector_mus_e_ho               , "museho"                );
  iEvent.put(vector_mus_e_emS9             , "museemS9"              );
  iEvent.put(vector_mus_e_hadS9            , "musehadS9"             );
  iEvent.put(vector_mus_e_hoS9             , "musehoS9"              );
  iEvent.put(vector_mus_iso_trckvetoDep    , "musisotrckvetoDep"     );
  iEvent.put(vector_mus_iso_ecalvetoDep    , "musisoecalvetoDep"     );
  iEvent.put(vector_mus_iso_hcalvetoDep    , "musisohcalvetoDep"     );
  iEvent.put(vector_mus_iso_hovetoDep      , "musisohovetoDep"       );
  iEvent.put(vector_mus_iso03_sumPt        , "musiso03sumPt"         );
  iEvent.put(vector_mus_iso03_emEt         , "musiso03emEt"          );
  iEvent.put(vector_mus_iso03_hadEt        , "musiso03hadEt"         );
  iEvent.put(vector_mus_iso03_hoEt         , "musiso03hoEt"          );
  iEvent.put(vector_mus_iso03_ntrk         , "musiso03ntrk"          );
  iEvent.put(vector_mus_iso05_sumPt        , "musiso05sumPt"         );
  iEvent.put(vector_mus_iso05_emEt         , "musiso05emEt"          );
  iEvent.put(vector_mus_iso05_hadEt        , "musiso05hadEt"         );
  iEvent.put(vector_mus_iso05_hoEt         , "musiso05hoEt"          );
  iEvent.put(vector_mus_iso05_ntrk         , "musiso05ntrk"          );

  
  iEvent.put(vector_mus_gfit_d0            , "musgfitd0"             );
  iEvent.put(vector_mus_gfit_z0            , "musgfitz0"             );
  iEvent.put(vector_mus_gfit_d0Err         , "musgfitd0Err"          );
  iEvent.put(vector_mus_gfit_z0Err         , "musgfitz0Err"          );
  iEvent.put(vector_mus_gfit_d0corr        , "musgfitd0corr"         );
  iEvent.put(vector_mus_gfit_z0corr        , "musgfitz0corr"         );
  iEvent.put(vector_mus_gfit_qoverp        , "musgfitqoverp"         );
  iEvent.put(vector_mus_gfit_qoverpError   , "musgfitqoverpError"    );
  iEvent.put(vector_mus_gfit_chi2          , "musgfitchi2"           );
  iEvent.put(vector_mus_gfit_ndof          , "musgfitndof"           );
  iEvent.put(vector_mus_gfit_validHits     , "musgfitvalidHits"      );

  iEvent.put(vector_mus_sta_d0             , "musstad0"              );
  iEvent.put(vector_mus_sta_z0             , "musstaz0"              );
  iEvent.put(vector_mus_sta_d0Err          , "musstad0Err"           );
  iEvent.put(vector_mus_sta_z0Err          , "musstaz0Err"           );
  iEvent.put(vector_mus_sta_d0corr         , "musstad0corr"          );
  iEvent.put(vector_mus_sta_z0corr         , "musstaz0corr"          );
  iEvent.put(vector_mus_sta_qoverp         , "musstaqoverp"          );
  iEvent.put(vector_mus_sta_qoverpError    , "musstaqoverpError"     );
  iEvent.put(vector_mus_sta_chi2           , "musstachi2"            );
  iEvent.put(vector_mus_sta_ndof           , "musstandof"            );
  iEvent.put(vector_mus_sta_validHits      , "musstavalidHits"       );

  
  iEvent.put(vector_mus_timeNumStationsUsed, "mustimeNumStationsUsed"); 
  iEvent.put(vector_mus_timeAtIpInOut      , "mustimeAtIpInOut"      );
  iEvent.put(vector_mus_timeAtIpInOutErr   , "mustimeAtIpInOutErr"   );
  iEvent.put(vector_mus_timeAtIpOutIn      , "mustimeAtIpOutIn"      );
  iEvent.put(vector_mus_timeAtIpOutInErr   , "mustimeAtIpOutInErr"   );
  iEvent.put(vector_mus_timeDirection      , "mustimeDirection"      );
  iEvent.put(vector_mus_pid_TMLastStationLoose	        , "muspidTMLastStationLoose"     );
  iEvent.put(vector_mus_pid_TMLastStationTight	        , "muspidTMLastStationTight"     );
  iEvent.put(vector_mus_pid_TM2DCompatibilityLoose	, "muspidTM2DCompatibilityLoose" );
  iEvent.put(vector_mus_pid_TM2DCompatibilityTight	, "muspidTM2DCompatibilityTight" );
  iEvent.put(vector_mus_caloCompatibility               , "muscaloCompatibility"         );
  iEvent.put(vector_mus_vertex_p4                       , "musvertexp4"                  );
  iEvent.put(vector_mus_gfit_vertex_p4                  , "musgfitvertexp4"              );
  iEvent.put(vector_mus_gfit_outerPos_p4                , "musgfitouterPosp4"            );
  iEvent.put(vector_mus_sta_vertex_p4                   , "musstavertexp4"               );
  iEvent.put(vector_mus_fitdefault_p4      , "musfitdefaultp4" );
  iEvent.put(vector_mus_fitfirsthit_p4     , "musfitfirsthitp4");
  iEvent.put(vector_mus_fitpicky_p4        , "musfitpickyp4"   );
  iEvent.put(vector_mus_fittev_p4          , "musfittevp4"     );

}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonMaker::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonMaker);
