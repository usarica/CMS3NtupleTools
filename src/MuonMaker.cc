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
// $Id: MuonMaker.cc,v 1.41 2010/05/06 20:10:08 kalavase Exp $
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

MuonMaker::MuonMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // mu track quantities
  produces<vector<int> >	     (branchprefix+"type"		      ).setBranchAlias(aliasprefix_+"_type"               );	// type
  produces<vector<int> >	     (branchprefix+"goodmask"	        ).setBranchAlias(aliasprefix_+"_goodmask"           ); // good mask
  produces<vector<LorentzVector> >   (branchprefix+"p4"		        ).setBranchAlias(aliasprefix_+"_p4"                 ); // candidate p4->this can either be gfit p4, tracker p4 or STA p4 (only for STA muoons)						
  produces<vector<LorentzVector> >   (branchprefix+"trkp4"	        ).setBranchAlias(aliasprefix_+"_trk_p4"             ); // track p4						
  produces<vector<LorentzVector> >   (branchprefix+"gfitp4"             ).setBranchAlias(aliasprefix_+"_gfit_p4"            ); // global fit p4, if global fit exists
  produces<vector<LorentzVector> >   (branchprefix+"stap4"              ).setBranchAlias(aliasprefix_+"_sta_p4"             ); // global fit p4, if global fit exists
  produces<vector<LorentzVector> >   (branchprefix+"ecalposp4"          ).setBranchAlias(aliasprefix_+"_ecalpos_p4"         ); // muon position at the ecal face
  produces<vector<int>   >           (branchprefix+"trkidx"             ).setBranchAlias(aliasprefix_+"_trkidx"             );	// track index matched to muon
  produces<vector<float> >	     (branchprefix+"d0"		      ).setBranchAlias(aliasprefix_+"_d0"                 ); // impact parameter at the point of closest approach	using the tracker fit
  produces<vector<float> >	     (branchprefix+"z0"		      ).setBranchAlias(aliasprefix_+"_z0"                 ); // z position of the point of closest approach. From the si track		
  produces<vector<float> >	     (branchprefix+"d0corr"	      ).setBranchAlias(aliasprefix_+"_d0corr"             ); // corrected impact parameter at the point of closest approach. From si track	
  produces<vector<float> >	     (branchprefix+"z0corr"	      ).setBranchAlias(aliasprefix_+"_z0corr"             ); // corrected z position of the point of closest approach. From si track		
  produces<vector<float> >	     (branchprefix+"vertexphi"	      ).setBranchAlias(aliasprefix_+"_vertexphi"          ); // phi angle of the point of closest approach. From si track		
  produces<vector<float> >	     (branchprefix+"chi2"		      ).setBranchAlias(aliasprefix_+"_chi2"               ); // chi2 of the silicon tracker fit			
  produces<vector<float> >	     (branchprefix+"ndof"		      ).setBranchAlias(aliasprefix_+"_ndof"               ); // number of degrees of freedom of the si tracker fit		
  produces<vector<int> >	     (branchprefix+"validHits"	      ).setBranchAlias(aliasprefix_+"_validHits"          ); // number of used hits in the sitracker fit			
  produces<vector<int> >	     (branchprefix+"lostHits"	      ).setBranchAlias(aliasprefix_+"_lostHits"           ); // number of lost hits in the sitracker fit			
  produces<vector<int> >             (branchprefix+"gfitvalidSTAHits"   ).setBranchAlias(aliasprefix_+"_gfit_validSTAHits"  ); // number of hits in the stand alone fit that made it into the gfit
  produces<vector<int> >             (branchprefix+"gfitvalidSiHits"    ).setBranchAlias(aliasprefix_+"_gfit_validSiHits"   ); // number of hits in the Si fit that made it into the gfit
  produces<vector<float> >	     (branchprefix+"d0Err"	      ).setBranchAlias(aliasprefix_+"_d0Err"              ); // error on the impact parameter, si track fit			
  produces<vector<float> >	     (branchprefix+"z0Err"	      ).setBranchAlias(aliasprefix_+"_z0Err"              ); // error on z position of the point of closest approach, si track fit	
  produces<vector<float> >	     (branchprefix+"ptErr"	      ).setBranchAlias(aliasprefix_+"_ptErr"              ); // si track Pt error					
  produces<vector<float> >	     (branchprefix+"etaErr"	      ).setBranchAlias(aliasprefix_+"_etaErr"             ); // si track eta error					
  produces<vector<float> >	     (branchprefix+"phiErr"	      ).setBranchAlias(aliasprefix_+"_phiErr"             ); // si track phi error					
  produces<vector<int> >	     (branchprefix+"charge"	      ).setBranchAlias(aliasprefix_+"_charge"             ); // charge from muon object 						
  produces<vector<int> >	     (branchprefix+"trkcharge"	      ).setBranchAlias(aliasprefix_+"_trk_charge"         ); // si track charge
  produces<vector<float> >           (branchprefix+"qoverp"             ).setBranchAlias(aliasprefix_+"_qoverp"             ); // si track qoverp
  produces<vector<float> >           (branchprefix+"qoverpError"        ).setBranchAlias(aliasprefix_+"_qoverpError"        ); // si track qoverp error
    // muon quantities
  produces<vector<int> >             (branchprefix+"nmatches"	      ).setBranchAlias(aliasprefix_+"_nmatches"           ); // number of stations with matched segments 
  produces<vector<float> >	     (branchprefix+"eem"		      ).setBranchAlias(aliasprefix_+"_e_em"               ); // energy in crossed ECAL crystalls 
  produces<vector<float> >	     (branchprefix+"ehad"		      ).setBranchAlias(aliasprefix_+"_e_had"              ); // energy in crossed HCAL towers 
  produces<vector<float> >	     (branchprefix+"eho"		      ).setBranchAlias(aliasprefix_+"_e_ho"               ); // energy in crossed HO towers 
  produces<vector<float> >	     (branchprefix+"eemS9"	      ).setBranchAlias(aliasprefix_+"_e_emS9"             ); // energy in 3x3 ECAL crystall shape 
  produces<vector<float> >	     (branchprefix+"ehadS9"	      ).setBranchAlias(aliasprefix_+"_e_hadS9"            ); //energy in 3x3 HCAL towers 
  produces<vector<float> >	     (branchprefix+"ehoS9"	      ).setBranchAlias(aliasprefix_+"_e_hoS9"             ); // energy in 3x3 HO towers 
  produces<vector<float> >           (branchprefix+"isotrckvetoDep"     ).setBranchAlias(aliasprefix_+"_iso_trckvetoDep"    );//sumPt in the veto cone, tracker
  produces<vector<float> >           (branchprefix+"isoecalvetoDep"     ).setBranchAlias(aliasprefix_+"_iso_ecalvetoDep"    );//sumEt in the veto cone, ecal
  produces<vector<float> >           (branchprefix+"isohcalvetoDep"     ).setBranchAlias(aliasprefix_+"_iso_hcalvetoDep"    );//sumPt in the veto cone, hcal
  produces<vector<float> >           (branchprefix+"isohovetoDep"       ).setBranchAlias(aliasprefix_+"_iso_hovetoDep"      );//sumPt in the veto cone, ho
  produces<vector<float> >	     (branchprefix+"iso03sumPt"	      ).setBranchAlias(aliasprefix_+"_iso03_sumPt"        ); // sum of track Pt for cone of 0.3 
  produces<vector<float> >	     (branchprefix+"iso03emEt"	      ).setBranchAlias(aliasprefix_+"_iso03_emEt"         ); // sum of ecal Et for cone of 0.3 
  produces<vector<float> >	     (branchprefix+"iso03hadEt"	      ).setBranchAlias(aliasprefix_+"_iso03_hadEt"        ); // sum of hcal Et for cone of 0.3 
  produces<vector<float> >	     (branchprefix+"iso03hoEt"	      ).setBranchAlias(aliasprefix_+"_iso03_hoEt"         ); // sum of ho Et for cone of 0.3 
  produces<vector<int> >	     (branchprefix+"iso03ntrk"	      ).setBranchAlias(aliasprefix_+"_iso03_ntrk"         ); // number of tracks in the cone of 0.3 
  produces<vector<float> >	     (branchprefix+"iso05sumPt"	      ).setBranchAlias(aliasprefix_+"_iso05_sumPt"        ); // sum of track Pt for cone of 0.5 
  produces<vector<float> >	     (branchprefix+"iso05emEt"          ).setBranchAlias(aliasprefix_+"_iso05_emEt"         ); // sum of ecal Et for cone of 0.5 
  produces<vector<float> >	     (branchprefix+"iso05hadEt"	      ).setBranchAlias(aliasprefix_+"_iso05_hadEt"        ); // sum of hcal Et for cone of 0.5 
  produces<vector<float> >	     (branchprefix+"iso05hoEt"	      ).setBranchAlias(aliasprefix_+"_iso05_hoEt"         ); // sum of ho Et for cone of 0.5 
  produces<vector<int> >	     (branchprefix+"iso05ntrk"	      ).setBranchAlias(aliasprefix_+"_iso05_ntrk"         ); // number of tracks in the cone of 0.5 

  //new
  produces<vector<float> >           (branchprefix+"gfitd0"             ).setBranchAlias(aliasprefix_+"_gfit_d0"            ); // d0 from global fit, if it exists
  produces<vector<float> >           (branchprefix+"gfitz0"             ).setBranchAlias(aliasprefix_+"_gfit_z0"            ); // z0 from global fit, if it exists
  produces<vector<float> >           (branchprefix+"gfitd0Err"          ).setBranchAlias(aliasprefix_+"_gfit_d0Err"         ); // d0Err from global fit, if it exists
  produces<vector<float> >           (branchprefix+"gfitz0Err"          ).setBranchAlias(aliasprefix_+"_gfit_z0Err"         ); // z0Err from global fit, if it exists
  produces<vector<float> >           (branchprefix+"gfitd0corr"         ).setBranchAlias(aliasprefix_+"_gfit_d0corr"        ); // Beamspot corrected d0 from global fit, if it exists
  produces<vector<float> >           (branchprefix+"gfitz0corr"         ).setBranchAlias(aliasprefix_+"_gfit_z0corr"        ); // Beamspot corrected z0 from global fit, if it exists
  produces<vector<float> >           (branchprefix+"gfitqoverp"         ).setBranchAlias(aliasprefix_+"_gfit_qoverp"        ); // global track qoverp
  produces<vector<float> >           (branchprefix+"gfitqoverpError"    ).setBranchAlias(aliasprefix_+"_gfit_qoverpError"   ); // global track qoverp error  
  
  produces<vector<float> >	     (branchprefix+"gfitchi2"           ).setBranchAlias(aliasprefix_+"_gfit_chi2"          ); // chi2 of the global muon fit 
  produces<vector<float> >	     (branchprefix+"gfitndof"	      ).setBranchAlias(aliasprefix_+"_gfit_ndof"          ); // number of degree of freedom of the global muon fit 
  produces<vector<int> >	     (branchprefix+"gfitvalidHits"      ).setBranchAlias(aliasprefix_+"_gfit_validHits"     ); // number of valid hits of the global muon fit 
  //STA fit crap
  produces<vector<float> >           (branchprefix+"stad0"             ).setBranchAlias(aliasprefix_+"_sta_d0"            ); // d0 from STA fit, if it exists
  produces<vector<float> >           (branchprefix+"staz0"             ).setBranchAlias(aliasprefix_+"_sta_z0"            ); // z0 from STA fit, if it exists
  produces<vector<float> >           (branchprefix+"stad0Err"          ).setBranchAlias(aliasprefix_+"_sta_d0Err"         ); // d0Err from STA fit, if it exists
  produces<vector<float> >           (branchprefix+"staz0Err"          ).setBranchAlias(aliasprefix_+"_sta_z0Err"         ); // z0Err from STA fit, if it exists
  produces<vector<float> >           (branchprefix+"stad0corr"         ).setBranchAlias(aliasprefix_+"_sta_d0corr"        ); // Beamspot corrected d0 from STA fit, if it exists
  produces<vector<float> >           (branchprefix+"staz0corr"         ).setBranchAlias(aliasprefix_+"_sta_z0corr"        ); // Beamspot corrected z0 from STA fit, if it exists
  produces<vector<float> >           (branchprefix+"staqoverp"         ).setBranchAlias(aliasprefix_+"_sta_qoverp"        ); // STA track qoverp
  produces<vector<float> >           (branchprefix+"staqoverpError"    ).setBranchAlias(aliasprefix_+"_sta_qoverpError"   ); // STA track qoverp error  
  
  produces<vector<float> >	     (branchprefix+"stachi2"           ).setBranchAlias(aliasprefix_+"_sta_chi2"          ); // chi2 of the STA muon fit 
  produces<vector<float> >	     (branchprefix+"standof"	     ).setBranchAlias(aliasprefix_+"_sta_ndof"          ); // number of degree of freedom of the STA muon fit 
  produces<vector<int> >	     (branchprefix+"stavalidHits"      ).setBranchAlias(aliasprefix_+"_sta_validHits"     ); // number of valid hits of the STA muon fit 
  

  
  //Muon timing info -> http://cmslxr.fnal.gov/lxr/source/DataFormats/MuonReco/interface/MuonTime.h
  produces<vector<int> >             (branchprefix+"timeNumStationsUsed").setBranchAlias(aliasprefix_+"_timeNumStationsUsed"); //number of muon stations used for timing info
  // time of arrival at the IP for the Beta=1 hypothesis -> particle moving from inside out
  produces<vector<float> >           (branchprefix+"timeAtIpInOut"      ).setBranchAlias(aliasprefix_+"_timeAtIpInOut"      ); 
  produces<vector<float> >           (branchprefix+"timeAtIpInOutErr"   ).setBranchAlias(aliasprefix_+"_timeAtIpInOutErr"   );
  //particle moving from outside in
  produces<vector<float> >           (branchprefix+"timeAtIpOutIn"      ).setBranchAlias(aliasprefix_+"_timeAtIpOutIn"      );
  produces<vector<float> >           (branchprefix+"timeAtIpOutInErr"   ).setBranchAlias(aliasprefix_+"_timeAtIpOutInErr"   );
  //direction estimate based on time dispersion. Enum defn given in above header file. in 312: 
  //Direction { OutsideIn = -1, Undefined = 0, InsideOut = 1 };
  produces<vector<int> >             (branchprefix+"timeDirection"      ).setBranchAlias(aliasprefix_+"_timeDirection"      );

  
  
  // loose tracker muon identification based on muon/hadron penetration depth difference       
  produces<vector<int> >	     (branchprefix+"pidTMLastStationLoose"    ).setBranchAlias(aliasprefix_+"_pid_TMLastStationLoose"    ); 
  // tight tracker muon identification based on muon/hadron penetration depth difference       
  produces<vector<int> >	     (branchprefix+"pidTMLastStationTight"    ).setBranchAlias(aliasprefix_+"_pid_TMLastStationTight"    );
  // loose tracker muon likelihood identification based on muon matches and calo depositions   
  produces<vector<int> >	     (branchprefix+"pidTM2DCompatibilityLoose").setBranchAlias(aliasprefix_+"_pid_TM2DCompatibilityLoose");
  // tight tracker muon likelihood identification based on muon matches and calo depositions
  produces<vector<int> >	     (branchprefix+"pidTM2DCompatibilityTight").setBranchAlias(aliasprefix_+"_pid_TM2DCompatibilityTight"); 
  //calo compatibility variable
  produces<vector<float> >           (branchprefix+"caloCompatibility").setBranchAlias(aliasprefix_+"_caloCompatibility");
  //overlap index (-1 if none)
  produces<vector<int> >           (branchprefix+"nOverlaps").setBranchAlias(aliasprefix_+"_nOverlaps");
  produces<vector<int> >           (branchprefix+"overlap0").setBranchAlias(aliasprefix_+"_overlap0");
  produces<vector<int> >           (branchprefix+"overlap1").setBranchAlias(aliasprefix_+"_overlap1");
  //p4 because we're not able to (yet) read XYZPointDs in bare root for some reason 
  //the 4th co-ordinate is 0
  produces<vector<LorentzVector> >   (branchprefix+"vertexp4"         ).setBranchAlias(aliasprefix_+"_vertex_p4" ); // from the silicon fit
  produces<vector<LorentzVector> >   (branchprefix+"gfitvertexp4"     ).setBranchAlias(aliasprefix_+"_gfit_vertex_p4");
  produces<vector<LorentzVector> >   (branchprefix+"gfitouterPosp4"   ).setBranchAlias(aliasprefix_+"_gfit_outerPos_p4");
  produces<vector<LorentzVector> >   (branchprefix+"stavertexp4"      ).setBranchAlias(aliasprefix_+"_sta_vertex_p4");
  produces<vector<LorentzVector> >   (branchprefix+"fitdefaultp4" ).setBranchAlias(aliasprefix_+"_fitdefault_p4" );
  produces<vector<LorentzVector> >   (branchprefix+"fitfirsthitp4").setBranchAlias(aliasprefix_+"_fitfirsthit_p4");
  produces<vector<LorentzVector> >   (branchprefix+"fitpickyp4"   ).setBranchAlias(aliasprefix_+"_fitpicky_p4"	 );
  produces<vector<LorentzVector> >   (branchprefix+"fittevp4"     ).setBranchAlias(aliasprefix_+"_fittev_p4"     );
  
  muonsInputTag  		= iConfig.getParameter<edm::InputTag>("muonsInputTag" ); 
  beamSpotInputTag 		= iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
  tevMuonsName    		= iConfig.getParameter<string>("tevMuonsName" ); 
 
}

void MuonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  // make vectors to hold the information
  auto_ptr<vector<int> >	   vector_mus_type    	          (new vector<int>	       );        
  auto_ptr<vector<int> >	   vector_mus_goodmask            (new vector<int>             );        
  auto_ptr<vector<LorentzVector> > vector_mus_p4                  (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > vector_mus_trk_p4	          (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > vector_mus_gfit_p4             (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > vector_mus_sta_p4              (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > vector_mus_ecalpos_p4          (new vector<LorentzVector>   );
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
  auto_ptr<vector<int> >           vector_mus_nOverlaps                    (new vector<int>  );
  auto_ptr<vector<int> >           vector_mus_overlap0                     (new vector<int>  );
  auto_ptr<vector<int> >           vector_mus_overlap1                     (new vector<int>  );

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
    
    for (int iG = 0; iG < 24; ++iG){ //overkill here
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
    math::XYZPoint ecal_p(-9999., -9999., -9999.);
    if(muon->isEnergyValid())
      ecal_p = muon->calEnergy().ecal_position;
    vector_mus_ecalpos_p4         ->push_back(LorentzVector(ecal_p.x(), ecal_p.y(), ecal_p.z(), 0.0)                               );
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

    int mus_overlap0 = -1;
    int mus_overlap1 = -1;
    int muInd = -1;
    int mus_nOverlaps = 0;
    for (edm::View<Muon>::const_iterator muonJ = muon_h->begin(); muonJ != muons_end; ++muonJ) {
      muInd++;
      if (muonJ!=muon){
	if (muon::overlap(*muon, *muonJ)){
	  if (mus_overlap0 == -1) mus_overlap0 = muInd;
	  if (mus_overlap0 != -1) mus_overlap1 = muInd;
	  mus_nOverlaps++;
	}
      }
    }
    
    vector_mus_nOverlaps          ->push_back(mus_nOverlaps );
    vector_mus_overlap0          ->push_back(mus_overlap0 );
    vector_mus_overlap1          ->push_back(mus_overlap1 );
    

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
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(vector_mus_type               , branchprefix+"type"               );
  iEvent.put(vector_mus_goodmask           , branchprefix+"goodmask"           );
  iEvent.put(vector_mus_p4                 , branchprefix+"p4"                 );
  iEvent.put(vector_mus_trk_p4             , branchprefix+"trkp4"              );
  iEvent.put(vector_mus_gfit_p4            , branchprefix+"gfitp4"             );
  iEvent.put(vector_mus_sta_p4             , branchprefix+"stap4"              );
  iEvent.put(vector_mus_ecalpos_p4         , branchprefix+"ecalposp4"          ); 
  iEvent.put(vector_mus_trkidx             , branchprefix+"trkidx"             );
  iEvent.put(vector_mus_d0                 , branchprefix+"d0"                 );
  iEvent.put(vector_mus_z0                 , branchprefix+"z0"                 );
  iEvent.put(vector_mus_d0corr             , branchprefix+"d0corr"             );
  iEvent.put(vector_mus_z0corr             , branchprefix+"z0corr"             );
  iEvent.put(vector_mus_vertexphi          , branchprefix+"vertexphi"          );
  iEvent.put(vector_mus_chi2               , branchprefix+"chi2"               );
  iEvent.put(vector_mus_ndof               , branchprefix+"ndof"               );
  iEvent.put(vector_mus_validHits          , branchprefix+"validHits"          );
  iEvent.put(vector_mus_lostHits           , branchprefix+"lostHits"           );
  iEvent.put(vector_mus_gfit_validSTAHits  , branchprefix+"gfitvalidSTAHits"   );
  iEvent.put(vector_mus_gfit_validSiHits   , branchprefix+"gfitvalidSiHits"    );
  iEvent.put(vector_mus_d0Err              , branchprefix+"d0Err"              );
  iEvent.put(vector_mus_z0Err              , branchprefix+"z0Err"              );
  iEvent.put(vector_mus_ptErr              , branchprefix+"ptErr"              );
  iEvent.put(vector_mus_etaErr             , branchprefix+"etaErr"             );
  iEvent.put(vector_mus_phiErr             , branchprefix+"phiErr"             );
  iEvent.put(vector_mus_charge             , branchprefix+"charge"             );
  iEvent.put(vector_mus_trk_charge         , branchprefix+"trkcharge"          );
  iEvent.put(vector_mus_qoverp             , branchprefix+"qoverp"             );
  iEvent.put(vector_mus_qoverpError        , branchprefix+"qoverpError"        );
  iEvent.put(vector_mus_nmatches           , branchprefix+"nmatches"           );
  iEvent.put(vector_mus_e_em               , branchprefix+"eem"                );
  iEvent.put(vector_mus_e_had              , branchprefix+"ehad"               );
  iEvent.put(vector_mus_e_ho               , branchprefix+"eho"                );
  iEvent.put(vector_mus_e_emS9             , branchprefix+"eemS9"              );
  iEvent.put(vector_mus_e_hadS9            , branchprefix+"ehadS9"             );
  iEvent.put(vector_mus_e_hoS9             , branchprefix+"ehoS9"              );
  iEvent.put(vector_mus_iso_trckvetoDep    , branchprefix+"isotrckvetoDep"     );
  iEvent.put(vector_mus_iso_ecalvetoDep    , branchprefix+"isoecalvetoDep"     );
  iEvent.put(vector_mus_iso_hcalvetoDep    , branchprefix+"isohcalvetoDep"     );
  iEvent.put(vector_mus_iso_hovetoDep      , branchprefix+"isohovetoDep"       );
  iEvent.put(vector_mus_iso03_sumPt        , branchprefix+"iso03sumPt"         );
  iEvent.put(vector_mus_iso03_emEt         , branchprefix+"iso03emEt"          );
  iEvent.put(vector_mus_iso03_hadEt        , branchprefix+"iso03hadEt"         );
  iEvent.put(vector_mus_iso03_hoEt         , branchprefix+"iso03hoEt"          );
  iEvent.put(vector_mus_iso03_ntrk         , branchprefix+"iso03ntrk"          );
  iEvent.put(vector_mus_iso05_sumPt        , branchprefix+"iso05sumPt"         );
  iEvent.put(vector_mus_iso05_emEt         , branchprefix+"iso05emEt"          );
  iEvent.put(vector_mus_iso05_hadEt        , branchprefix+"iso05hadEt"         );
  iEvent.put(vector_mus_iso05_hoEt         , branchprefix+"iso05hoEt"          );
  iEvent.put(vector_mus_iso05_ntrk         , branchprefix+"iso05ntrk"          );

  
  iEvent.put(vector_mus_gfit_d0            , branchprefix+"gfitd0"             );
  iEvent.put(vector_mus_gfit_z0            , branchprefix+"gfitz0"             );
  iEvent.put(vector_mus_gfit_d0Err         , branchprefix+"gfitd0Err"          );
  iEvent.put(vector_mus_gfit_z0Err         , branchprefix+"gfitz0Err"          );
  iEvent.put(vector_mus_gfit_d0corr        , branchprefix+"gfitd0corr"         );
  iEvent.put(vector_mus_gfit_z0corr        , branchprefix+"gfitz0corr"         );
  iEvent.put(vector_mus_gfit_qoverp        , branchprefix+"gfitqoverp"         );
  iEvent.put(vector_mus_gfit_qoverpError   , branchprefix+"gfitqoverpError"    );
  iEvent.put(vector_mus_gfit_chi2          , branchprefix+"gfitchi2"           );
  iEvent.put(vector_mus_gfit_ndof          , branchprefix+"gfitndof"           );
  iEvent.put(vector_mus_gfit_validHits     , branchprefix+"gfitvalidHits"      );

  iEvent.put(vector_mus_sta_d0             , branchprefix+"stad0"              );
  iEvent.put(vector_mus_sta_z0             , branchprefix+"staz0"              );
  iEvent.put(vector_mus_sta_d0Err          , branchprefix+"stad0Err"           );
  iEvent.put(vector_mus_sta_z0Err          , branchprefix+"staz0Err"           );
  iEvent.put(vector_mus_sta_d0corr         , branchprefix+"stad0corr"          );
  iEvent.put(vector_mus_sta_z0corr         , branchprefix+"staz0corr"          );
  iEvent.put(vector_mus_sta_qoverp         , branchprefix+"staqoverp"          );
  iEvent.put(vector_mus_sta_qoverpError    , branchprefix+"staqoverpError"     );
  iEvent.put(vector_mus_sta_chi2           , branchprefix+"stachi2"            );
  iEvent.put(vector_mus_sta_ndof           , branchprefix+"standof"            );
  iEvent.put(vector_mus_sta_validHits      , branchprefix+"stavalidHits"       );

  
  iEvent.put(vector_mus_timeNumStationsUsed, branchprefix+"timeNumStationsUsed"); 
  iEvent.put(vector_mus_timeAtIpInOut      , branchprefix+"timeAtIpInOut"      );
  iEvent.put(vector_mus_timeAtIpInOutErr   , branchprefix+"timeAtIpInOutErr"   );
  iEvent.put(vector_mus_timeAtIpOutIn      , branchprefix+"timeAtIpOutIn"      );
  iEvent.put(vector_mus_timeAtIpOutInErr   , branchprefix+"timeAtIpOutInErr"   );
  iEvent.put(vector_mus_timeDirection      , branchprefix+"timeDirection"      );
  iEvent.put(vector_mus_pid_TMLastStationLoose	        , branchprefix+"pidTMLastStationLoose"     );
  iEvent.put(vector_mus_pid_TMLastStationTight	        , branchprefix+"pidTMLastStationTight"     );
  iEvent.put(vector_mus_pid_TM2DCompatibilityLoose	, branchprefix+"pidTM2DCompatibilityLoose" );
  iEvent.put(vector_mus_pid_TM2DCompatibilityTight	, branchprefix+"pidTM2DCompatibilityTight" );
  iEvent.put(vector_mus_caloCompatibility               , branchprefix+"caloCompatibility"         );
  iEvent.put(vector_mus_nOverlaps                       , branchprefix+"nOverlaps"                 );
  iEvent.put(vector_mus_overlap0                        , branchprefix+"overlap0"                 );
  iEvent.put(vector_mus_overlap1                        , branchprefix+"overlap1"                 );

  iEvent.put(vector_mus_vertex_p4                       , branchprefix+"vertexp4"                  );
  iEvent.put(vector_mus_gfit_vertex_p4                  , branchprefix+"gfitvertexp4"              );
  iEvent.put(vector_mus_gfit_outerPos_p4                , branchprefix+"gfitouterPosp4"            );
  iEvent.put(vector_mus_sta_vertex_p4                   , branchprefix+"stavertexp4"               );
  iEvent.put(vector_mus_fitdefault_p4      , branchprefix+"fitdefaultp4" );
  iEvent.put(vector_mus_fitfirsthit_p4     , branchprefix+"fitfirsthitp4");
  iEvent.put(vector_mus_fitpicky_p4        , branchprefix+"fitpickyp4"   );
  iEvent.put(vector_mus_fittev_p4          , branchprefix+"fittevp4"     );

}


// ------------ method called once each job just before starting event loop  ------------
void MuonMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void MuonMaker::endJob() {}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonMaker);
