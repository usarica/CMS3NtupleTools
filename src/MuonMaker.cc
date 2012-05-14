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
// $Id: MuonMaker.cc,v 1.66 2012/05/14 17:15:16 cerati Exp $
//
//


// system include files
#include <memory>
#include <sstream>

// user include files
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/VectorUtil.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CMS2/NtupleMaker/interface/MuonMaker.h"


//////////////
// typedefs //
//////////////

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;


////////////////
// namespaces //
////////////////

using namespace std;
using namespace reco;
using namespace edm;


/////////////////
// Constructor //
/////////////////

MuonMaker::MuonMaker( const ParameterSet& iConfig ) {

  /////////////////////////////
  // Branch & Alias prefixes //
  /////////////////////////////

  aliasprefix_        = iConfig.getUntrackedParameter<string>("aliasPrefix");
  branchprefix_       = aliasprefix_; if( branchprefix_.find("_") != string::npos ) branchprefix_.replace( branchprefix_.find("_"), 1, "" );


  //////////////////////
  // Input Parameters //
  //////////////////////

  muonsInputTag    = iConfig.getParameter<InputTag> ("muonsInputTag"   );
  beamSpotInputTag = iConfig.getParameter<InputTag> ("beamSpotInputTag");
  pfCandsInputTag  = iConfig.getParameter<InputTag> ("pfCandsInputTag" );
  vtxInputTag      = iConfig.getParameter<InputTag> ("vtxInputTag"     );
  tevMuonsName     = iConfig.getParameter<string>   ("tevMuonsName"    );
  src_             = iConfig.getParameter<InputTag> ("cosmicCompat"    ); 


  /////////
  // STA //
  ///////// 

  produces<vector<float> >          ( branchprefix_ + "stad0"                     ).setBranchAlias( aliasprefix_ + "_sta_d0"            ); // d0 from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "staz0"                     ).setBranchAlias( aliasprefix_ + "_sta_z0"            ); // z0 from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "stad0Err"                  ).setBranchAlias( aliasprefix_ + "_sta_d0Err"         ); // d0Err from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "staz0Err"                  ).setBranchAlias( aliasprefix_ + "_sta_z0Err"         ); // z0Err from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "stad0corr"                 ).setBranchAlias( aliasprefix_ + "_sta_d0corr"        ); // Beamspot corrected d0 from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "staz0corr"                 ).setBranchAlias( aliasprefix_ + "_sta_z0corr"        ); // Beamspot corrected z0 from STA fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "staqoverp"                 ).setBranchAlias( aliasprefix_ + "_sta_qoverp"        ); // STA track qoverp
  produces<vector<float> >          ( branchprefix_ + "staqoverpError"            ).setBranchAlias( aliasprefix_ + "_sta_qoverpError"   ); // STA track qoverp error  
  produces<vector<float> >          ( branchprefix_ + "stachi2"                   ).setBranchAlias( aliasprefix_ + "_sta_chi2"          ); // chi2 of the STA muon fit 
  produces<vector<float> >          ( branchprefix_ + "standof"                   ).setBranchAlias( aliasprefix_ + "_sta_ndof"          ); // number of degree of freedom of the STA muon fit 
  produces<vector<int> >            ( branchprefix_ + "stavalidHits"              ).setBranchAlias( aliasprefix_ + "_sta_validHits"     ); // number of valid hits of the STA muon fit 
  produces<vector<LorentzVector> >  ( branchprefix_ + "stap4"                     ).setBranchAlias( aliasprefix_ + "_sta_p4"            ); // 
  produces<vector<LorentzVector> >  ( branchprefix_ + "stavertexp4"               ).setBranchAlias( aliasprefix_ + "_sta_vertex_p4"     );


  ////////////
  // Global //
  ////////////

  produces<vector<float> >          ( branchprefix_ + "gfitd0"                    ).setBranchAlias( aliasprefix_ + "_gfit_d0"            ); // d0 from global fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "gfitz0"                    ).setBranchAlias( aliasprefix_ + "_gfit_z0"            ); // z0 from global fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "gfitd0Err"                 ).setBranchAlias( aliasprefix_ + "_gfit_d0Err"         ); // d0Err from global fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "gfitz0Err"                 ).setBranchAlias( aliasprefix_ + "_gfit_z0Err"         ); // z0Err from global fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "gfitd0corr"                ).setBranchAlias( aliasprefix_ + "_gfit_d0corr"        ); // Beamspot corrected d0 from global fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "gfitz0corr"                ).setBranchAlias( aliasprefix_ + "_gfit_z0corr"        ); // Beamspot corrected z0 from global fit, if it exists
  produces<vector<float> >          ( branchprefix_ + "gfitqoverp"                ).setBranchAlias( aliasprefix_ + "_gfit_qoverp"        ); // global track qoverp
  produces<vector<float> >          ( branchprefix_ + "gfitqoverpError"           ).setBranchAlias( aliasprefix_ + "_gfit_qoverpError"   ); // global track qoverp error  
  produces<vector<float> >          ( branchprefix_ + "gfitchi2"                  ).setBranchAlias( aliasprefix_ + "_gfit_chi2"          ); // chi2 of the global muon fit 
  produces<vector<float> >          ( branchprefix_ + "gfitndof"                  ).setBranchAlias( aliasprefix_ + "_gfit_ndof"          ); // number of degree of freedom of the global muon fit 
  produces<vector<int> >            ( branchprefix_ + "gfitvalidHits"             ).setBranchAlias( aliasprefix_ + "_gfit_validHits"     ); // number of valid hits of the global muon fit 
  produces<vector<int> >            ( branchprefix_ + "gfitvalidSTAHits"          ).setBranchAlias( aliasprefix_ + "_gfit_validSTAHits"  ); // number of hits in the stand alone fit that made it into the gfit
  produces<vector<int> >            ( branchprefix_ + "gfitvalidSiHits"           ).setBranchAlias( aliasprefix_ + "_gfit_validSiHits"   ); // number of hits in the Si fit that made it into the gfit
  produces<vector<LorentzVector> >  ( branchprefix_ + "gfitp4"                    ).setBranchAlias( aliasprefix_ + "_gfit_p4"            ); // global fit p4, if global fit exists
  produces<vector<LorentzVector> >  ( branchprefix_ + "gfitvertexp4"              ).setBranchAlias( aliasprefix_ + "_gfit_vertex_p4"     );
  produces<vector<LorentzVector> >  ( branchprefix_ + "gfitouterPosp4"            ).setBranchAlias( aliasprefix_ + "_gfit_outerPos_p4"   );


  /////////////
  // Quality //
  /////////////

  produces<vector<bool> >           ( branchprefix_ + "updatedSta"                ).setBranchAlias( aliasprefix_ + "_updatedSta"          );  // Muon Quality - updatedSta
  produces<vector<bool> >           ( branchprefix_ + "tightMatch"                ).setBranchAlias( aliasprefix_ + "_tightMatch"          );  // Muon Quality - tightMatch
  produces<vector<float> >          ( branchprefix_ + "trkKink"                   ).setBranchAlias( aliasprefix_ + "_trkKink"             );  // Muon Quality - trkKink
  produces<vector<float> >          ( branchprefix_ + "glbKink"                   ).setBranchAlias( aliasprefix_ + "_glbKink"             );  // Muon Quality - glbKink
  produces<vector<float> >          ( branchprefix_ + "trkRelChi2"                ).setBranchAlias( aliasprefix_ + "_trkRelChi2"          );  // Muon Quality - trkRelChi2
  produces<vector<float> >          ( branchprefix_ + "staRelChi2"                ).setBranchAlias( aliasprefix_ + "_staRelChi2"          );  // Muon Quality - staRelChi2
  produces<vector<float> >          ( branchprefix_ + "chi2LocalPosition"         ).setBranchAlias( aliasprefix_ + "_chi2LocalPosition"   );  // Muon Quality - chi2LocalPositions
  produces<vector<float> >          ( branchprefix_ + "chi2LocalMomentum"         ).setBranchAlias( aliasprefix_ + "_chi2LocalMomentum"   );  // Muon Quality - chi2LocalMomentum
  produces<vector<float> >          ( branchprefix_ + "localDistance"             ).setBranchAlias( aliasprefix_ + "_localDistance"       );  // Muon Quality - localDistance
  produces<vector<float> >          ( branchprefix_ + "globalDeltaEtaPhi"         ).setBranchAlias( aliasprefix_ + "_globalDeltaEtaPhi"   );  // Muon Quality - globalDeltaEtaPhi
  produces<vector<float> >          ( branchprefix_ + "glbTrackProbability"       ).setBranchAlias( aliasprefix_ + "_glbTrackProbability" );  // Muon Quality - glbTrackProbability


  ////////////
  // Timing //
  //////////// 

  produces<vector<int> >            ( branchprefix_ + "timeNumStationsUsed"       ).setBranchAlias( aliasprefix_ + "_timeNumStationsUsed" ); // number of muon stations used for timing info
  produces<vector<float> >          ( branchprefix_ + "timeAtIpInOut"             ).setBranchAlias( aliasprefix_ + "_timeAtIpInOut"       ); // time of arrival at the IP for the Beta=1 hypothesis -> particle moving from inside out
  produces<vector<float> >          ( branchprefix_ + "timeAtIpInOutErr"          ).setBranchAlias( aliasprefix_ + "_timeAtIpInOutErr"    ); // particle moving from outside in
  produces<vector<float> >          ( branchprefix_ + "timeAtIpOutIn"             ).setBranchAlias( aliasprefix_ + "_timeAtIpOutIn"       );
  produces<vector<float> >          ( branchprefix_ + "timeAtIpOutInErr"          ).setBranchAlias( aliasprefix_ + "_timeAtIpOutInErr"    );
  produces<vector<int> >            ( branchprefix_ + "timeDirection"             ).setBranchAlias( aliasprefix_ + "_timeDirection"       ); //Direction { OutsideIn = -1, Undefined = 0, InsideOut = 1 };

                                    
  //////////////////////////
  // Cosmic Compatibility //
  //////////////////////////

  produces<vector<float> >          ( branchprefix_ + "cosmicCompat"              ).setBranchAlias( aliasprefix_ + "_cosmicCompat"        );
  produces<vector<float> >          ( branchprefix_ + "timeCompat"                ).setBranchAlias( aliasprefix_ + "_timeCompat"          );
  produces<vector<float> >          ( branchprefix_ + "backToBackCompat"          ).setBranchAlias( aliasprefix_ + "_backToBackCompat"    );
  produces<vector<float> >          ( branchprefix_ + "overlapCompat"             ).setBranchAlias( aliasprefix_ + "_overlapCompat"       );
  produces<vector<float> >          ( branchprefix_ + "vertexCompat"              ).setBranchAlias( aliasprefix_ + "_vertexCompat"        );


  //////////
  // Muon //
  //////////

  produces<vector<int> >            ( branchprefix_ + "type"                      ).setBranchAlias( aliasprefix_ + "_type"                   ); // type
  produces<vector<int> >            ( branchprefix_ + "goodmask"                  ).setBranchAlias( aliasprefix_ + "_goodmask"               ); // good mask
  produces<vector<int> >            ( branchprefix_ + "charge"                    ).setBranchAlias( aliasprefix_ + "_charge"                 ); // charge from muon object             
  produces<vector<int> >            ( branchprefix_ + "nmatches"                  ).setBranchAlias( aliasprefix_ + "_nmatches"               ); // number of stations with matched segments 
  produces<vector<int> >            ( branchprefix_ + "nOverlaps"                 ).setBranchAlias( aliasprefix_ + "_nOverlaps"              ); //overlap index (-1 if none)
  produces<vector<int> >            ( branchprefix_ + "overlap0"                  ).setBranchAlias( aliasprefix_ + "_overlap0"               );
  produces<vector<int> >            ( branchprefix_ + "overlap1"                  ).setBranchAlias( aliasprefix_ + "_overlap1"               );
  produces<vector<float> >          ( branchprefix_ + "caloCompatibility"         ).setBranchAlias( aliasprefix_ + "_caloCompatibility"      ); // calo compatibility variable
  produces<vector<float> >          ( branchprefix_ + "segmCompatibility"         ).setBranchAlias( aliasprefix_ + "_segmCompatibility"      );
  produces<vector<LorentzVector> >  ( branchprefix_ + "p4"                        ).setBranchAlias( aliasprefix_ + "_p4"                     ); // candidate p4->this can either be gfit p4, tracker p4 or STA p4 (only for STA muoons)     
  produces<vector<LorentzVector> >  ( branchprefix_ + "ecalposp4"                 ).setBranchAlias( aliasprefix_ + "_ecalpos_p4"             ); // muon position at the ecal face
  produces<vector<int> >            ( branchprefix_ + "numberOfMatchedStations"   ).setBranchAlias( aliasprefix_ + "_numberOfMatchedStations"); // number of muon stations with muon segements used in the fit

  ////////
  // ID //
  ////////

  produces<vector<int> >            ( branchprefix_ + "pidTMLastStationLoose"     ).setBranchAlias( aliasprefix_ + "_pid_TMLastStationLoose"     ); // loose tracker muon identification based on muon/hadron penetration depth difference       
  produces<vector<int> >            ( branchprefix_ + "pidTMLastStationTight"     ).setBranchAlias( aliasprefix_ + "_pid_TMLastStationTight"     ); // tight tracker muon identification based on muon/hadron penetration depth difference       
  produces<vector<int> >            ( branchprefix_ + "pidTM2DCompatibilityLoose" ).setBranchAlias( aliasprefix_ + "_pid_TM2DCompatibilityLoose" ); // loose tracker muon likelihood identification based on muon matches and calo depositions   
  produces<vector<int> >            ( branchprefix_ + "pidTM2DCompatibilityTight" ).setBranchAlias( aliasprefix_ + "_pid_TM2DCompatibilityTight" ); // tight tracker muon likelihood identification based on muon matches and calo depositions
  produces<vector<int> >            ( branchprefix_ + "pidPFMuon" ).setBranchAlias( aliasprefix_ + "_pid_PFMuon" ); // is particle flow muon

  ////////////
  // Energy //
  ////////////

  produces<vector<float> >          ( branchprefix_ + "eem"                       ).setBranchAlias( aliasprefix_ + "_e_em"                ); // energy in crossed ECAL crystalls 
  produces<vector<float> >          ( branchprefix_ + "ehad"                      ).setBranchAlias( aliasprefix_ + "_e_had"               ); // energy in crossed HCAL towers 
  produces<vector<float> >          ( branchprefix_ + "eho"                       ).setBranchAlias( aliasprefix_ + "_e_ho"                ); // energy in crossed HO towers 
  produces<vector<float> >          ( branchprefix_ + "eemS9"                     ).setBranchAlias( aliasprefix_ + "_e_emS9"              ); // energy in 3x3 ECAL crystall shape 
  produces<vector<float> >          ( branchprefix_ + "ehadS9"                    ).setBranchAlias( aliasprefix_ + "_e_hadS9"             ); // energy in 3x3 HCAL towers 
  produces<vector<float> >          ( branchprefix_ + "ehoS9"                     ).setBranchAlias( aliasprefix_ + "_e_hoS9"              ); // energy in 3x3 HO towers 


  ///////////////
  // Isolation //
  ///////////////
                                    
  produces<vector<float> >          ( branchprefix_ + "isotrckvetoDep"            ).setBranchAlias( aliasprefix_ + "_iso_trckvetoDep"     ); // sumPt in the veto cone, tracker
  produces<vector<float> >          ( branchprefix_ + "isoecalvetoDep"            ).setBranchAlias( aliasprefix_ + "_iso_ecalvetoDep"     ); // sumEt in the veto cone, ecal
  produces<vector<float> >          ( branchprefix_ + "isohcalvetoDep"            ).setBranchAlias( aliasprefix_ + "_iso_hcalvetoDep"     ); // sumPt in the veto cone, hcal
  produces<vector<float> >          ( branchprefix_ + "isohovetoDep"              ).setBranchAlias( aliasprefix_ + "_iso_hovetoDep"       ); // sumPt in the veto cone, ho
  produces<vector<float> >          ( branchprefix_ + "iso03sumPt"                ).setBranchAlias( aliasprefix_ + "_iso03_sumPt"         ); // sum of track Pt for cone of 0.3 
  produces<vector<float> >          ( branchprefix_ + "iso03emEt"                 ).setBranchAlias( aliasprefix_ + "_iso03_emEt"          ); // sum of ecal Et for cone of 0.3 
  produces<vector<float> >          ( branchprefix_ + "iso03hadEt"                ).setBranchAlias( aliasprefix_ + "_iso03_hadEt"         ); // sum of hcal Et for cone of 0.3 
  produces<vector<float> >          ( branchprefix_ + "iso03hoEt"                 ).setBranchAlias( aliasprefix_ + "_iso03_hoEt"          ); // sum of ho Et for cone of 0.3 
  produces<vector<int> >            ( branchprefix_ + "iso03ntrk"                 ).setBranchAlias( aliasprefix_ + "_iso03_ntrk"          ); // number of tracks in the cone of 0.3 
  produces<vector<float> >          ( branchprefix_ + "iso05sumPt"                ).setBranchAlias( aliasprefix_ + "_iso05_sumPt"         ); // sum of track Pt for cone of 0.5 
  produces<vector<float> >          ( branchprefix_ + "iso05emEt"                 ).setBranchAlias( aliasprefix_ + "_iso05_emEt"          ); // sum of ecal Et for cone of 0.5 
  produces<vector<float> >          ( branchprefix_ + "iso05hadEt"                ).setBranchAlias( aliasprefix_ + "_iso05_hadEt"         ); // sum of hcal Et for cone of 0.5 
  produces<vector<float> >          ( branchprefix_ + "iso05hoEt"                 ).setBranchAlias( aliasprefix_ + "_iso05_hoEt"          ); // sum of ho Et for cone of 0.5 
  produces<vector<int> >            ( branchprefix_ + "iso05ntrk"                 ).setBranchAlias( aliasprefix_ + "_iso05_ntrk"          ); // number of tracks in the cone of 0.5 


  ////////////
  // Tracks //
  ////////////

  produces<vector<int> >            ( branchprefix_ + "muonBestTrackType"         ).setBranchAlias( aliasprefix_ + "_muonBestTrackType"   ); 
  produces<vector<LorentzVector> >  ( branchprefix_ + "trkp4"                     ).setBranchAlias( aliasprefix_ + "_trk_p4"              ); // track p4            
  produces<vector<LorentzVector> >  ( branchprefix_ + "vertexp4"                  ).setBranchAlias( aliasprefix_ + "_vertex_p4"           ); // from the silicon fit
  produces<vector<int>   >          ( branchprefix_ + "trkidx"                    ).setBranchAlias( aliasprefix_ + "_trkidx"              ); // track index matched to muon
  produces<vector<float> >          ( branchprefix_ + "d0"                        ).setBranchAlias( aliasprefix_ + "_d0"                  ); // impact parameter at the point of closest approach  using the tracker fit
  produces<vector<float> >          ( branchprefix_ + "z0"                        ).setBranchAlias( aliasprefix_ + "_z0"                  ); // z position of the point of closest approach. From the si track    
  produces<vector<float> >          ( branchprefix_ + "d0corr"                    ).setBranchAlias( aliasprefix_ + "_d0corr"              ); // corrected impact parameter at the point of closest approach. From si track  
  produces<vector<float> >          ( branchprefix_ + "z0corr"                    ).setBranchAlias( aliasprefix_ + "_z0corr"              ); // corrected z position of the point of closest approach. From si track    
  produces<vector<float> >          ( branchprefix_ + "vertexphi"                 ).setBranchAlias( aliasprefix_ + "_vertexphi"           ); // phi angle of the point of closest approach. From si track    
  produces<vector<float> >          ( branchprefix_ + "chi2"                      ).setBranchAlias( aliasprefix_ + "_chi2"                ); // chi2 of the silicon tracker fit      
  produces<vector<float> >          ( branchprefix_ + "ndof"                      ).setBranchAlias( aliasprefix_ + "_ndof"                ); // number of degrees of freedom of the si tracker fit    
  produces<vector<int> >            ( branchprefix_ + "validHits"                 ).setBranchAlias( aliasprefix_ + "_validHits"           ); // number of used hits in the sitracker fit      
  produces<vector<int> >            ( branchprefix_ + "lostHits"                  ).setBranchAlias( aliasprefix_ + "_lostHits"            ); // number of lost hits in the sitracker fit      
  produces<vector<float> >          ( branchprefix_ + "d0Err"                     ).setBranchAlias( aliasprefix_ + "_d0Err"               ); // error on the impact parameter, si track fit      
  produces<vector<float> >          ( branchprefix_ + "z0Err"                     ).setBranchAlias( aliasprefix_ + "_z0Err"               ); // error on z position of the point of closest approach, si track fit  
  produces<vector<float> >          ( branchprefix_ + "ptErr"                     ).setBranchAlias( aliasprefix_ + "_ptErr"               ); // si track Pt error          
  produces<vector<float> >          ( branchprefix_ + "etaErr"                    ).setBranchAlias( aliasprefix_ + "_etaErr"              ); // si track eta error          
  produces<vector<float> >          ( branchprefix_ + "phiErr"                    ).setBranchAlias( aliasprefix_ + "_phiErr"              ); // si track phi error          
  produces<vector<int> >            ( branchprefix_ + "trkcharge"                 ).setBranchAlias( aliasprefix_ + "_trk_charge"          ); // si track charge
  produces<vector<float> >          ( branchprefix_ + "qoverp"                    ).setBranchAlias( aliasprefix_ + "_qoverp"              ); // si track qoverp
  produces<vector<float> >          ( branchprefix_ + "qoverpError"               ).setBranchAlias( aliasprefix_ + "_qoverpError"         ); // si track qoverp error

  produces<vector<LorentzVector> >  ( branchprefix_ + "fittpfmsp4"                ).setBranchAlias( aliasprefix_ + "_fittpfms_p4"         );
  produces<vector<LorentzVector> >  ( branchprefix_ + "fitpickyp4"                ).setBranchAlias( aliasprefix_ + "_fitpicky_p4"         );
  produces<vector<LorentzVector> >  ( branchprefix_ + "fittevp4"                  ).setBranchAlias( aliasprefix_ + "_fittev_p4"           );

  
  ////////
  // PF //
  ////////

  produces<vector< int> >           ( branchprefix_ + "pfcharge"                  ).setBranchAlias( aliasprefix_ + "_pfcharge"            );
  produces<vector< int> >           ( branchprefix_ + "pfparticleId"              ).setBranchAlias( aliasprefix_ + "_pfparticleId"        );
  produces<vector< int> >           ( branchprefix_ + "pfflag"                    ).setBranchAlias( aliasprefix_ + "_pfflag"              );
  produces<vector< float> >         ( branchprefix_ + "pfecalE"                   ).setBranchAlias( aliasprefix_ + "_pfecalE"             );
  produces<vector< float> >         ( branchprefix_ + "pfhcalE"                   ).setBranchAlias( aliasprefix_ + "_pfhcalE"             );
  produces<vector< float> >         ( branchprefix_ + "pfrawEcalE"                ).setBranchAlias( aliasprefix_ + "_pfrawEcalE"          );
  produces<vector< float> >         ( branchprefix_ + "pfrawHcalE"                ).setBranchAlias( aliasprefix_ + "_pfrawHcalE"          );
  produces<vector< float> >         ( branchprefix_ + "pfpS1E"                    ).setBranchAlias( aliasprefix_ + "_pfpS1E"              );
  produces<vector< float> >         ( branchprefix_ + "pfpS2E"                    ).setBranchAlias( aliasprefix_ + "_pfpS2E"              );
  produces<vector< float> >         ( branchprefix_ + "pfdeltaP"                  ).setBranchAlias( aliasprefix_ + "_pfdeltaP"            );
  produces<vector< float> >         ( branchprefix_ + "pfmvaepi"                  ).setBranchAlias( aliasprefix_ + "_pfmva_epi"           );
  produces<vector< float> >         ( branchprefix_ + "pfmvaemu"                  ).setBranchAlias( aliasprefix_ + "_pfmva_emu"           );
  produces<vector< float> >         ( branchprefix_ + "pfmvapimu"                 ).setBranchAlias( aliasprefix_ + "_pfmva_pimu"          );
  produces<vector< float> >         ( branchprefix_ + "pfmvanothinggamma"         ).setBranchAlias( aliasprefix_ + "_pfmva_nothing_gamma" );
  produces<vector< float> >         ( branchprefix_ + "pfmvanothingnh"            ).setBranchAlias( aliasprefix_ + "_pfmva_nothing_nh"    );
  produces<vector< LorentzVector> > ( branchprefix_ + "pfp4"                      ).setBranchAlias( aliasprefix_ + "_pfp4"                );
  produces<vector< LorentzVector> > ( branchprefix_ + "pfposAtEcalp4"             ).setBranchAlias( aliasprefix_ + "_pfposAtEcal_p4"      );
  /*
  produces<vector< float> >         ( branchprefix_ + "pfisoChargedHadrons"       ).setBranchAlias( aliasprefix_ + "_pfiso03ChargedHadrons");
  produces<vector< float> >         ( branchprefix_ + "pfisoNeutralHadrons"       ).setBranchAlias( aliasprefix_ + "_pfiso03NeutralHadrons");
  produces<vector< float> >         ( branchprefix_ + "pfisoPhotons"              ).setBranchAlias( aliasprefix_ + "_pfiso03Photons"       );
  produces<vector< float> >         ( branchprefix_ + "pfiso04ChargedHadrons"     ).setBranchAlias( aliasprefix_ + "_pfiso04ChargedHadrons");
  produces<vector< float> >         ( branchprefix_ + "pfiso04NeutralHadrons"     ).setBranchAlias( aliasprefix_ + "_pfiso04NeutralHadrons");
  produces<vector< float> >         ( branchprefix_ + "pfiso04Photons"            ).setBranchAlias( aliasprefix_ + "_pfiso04Photons"       );
  */
  produces<vector<float> >          ( branchprefix_ + "isoR03pfChargedHadronPt"   ).setBranchAlias( aliasprefix_ + "_isoR03_pf_ChargedHadronPt"   );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfChargedParticlePt" ).setBranchAlias( aliasprefix_ + "_isoR03_pf_ChargedParticlePt" );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfNeutralHadronEt"   ).setBranchAlias( aliasprefix_ + "_isoR03_pf_NeutralHadronEt"   );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfPhotonEt"          ).setBranchAlias( aliasprefix_ + "_isoR03_pf_PhotonEt"          );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfNeutralHadronEtHighThreshold").setBranchAlias( aliasprefix_ + "_isoR03_pf_NeutralHadronEtHighThreshold");
  produces<vector<float> >          ( branchprefix_ + "isoR03pfPhotonEtHighThreshold"       ).setBranchAlias( aliasprefix_ + "_isoR03_pf_PhotonEtHighThreshold"       );
  produces<vector<float> >          ( branchprefix_ + "isoR03pfPUPt"              ).setBranchAlias( aliasprefix_ + "_isoR03_pf_PUPt"              );

  produces<vector<float> >          ( branchprefix_ + "isoR04pfChargedHadronPt"   ).setBranchAlias( aliasprefix_ + "_isoR04_pf_ChargedHadronPt"   );
  produces<vector<float> >          ( branchprefix_ + "isoR04pfChargedParticlePt" ).setBranchAlias( aliasprefix_ + "_isoR04_pf_ChargedParticlePt" );
  produces<vector<float> >          ( branchprefix_ + "isoR04pfNeutralHadronEt"   ).setBranchAlias( aliasprefix_ + "_isoR04_pf_NeutralHadronEt"   );
  produces<vector<float> >          ( branchprefix_ + "isoR04pfPhotonEt"          ).setBranchAlias( aliasprefix_ + "_isoR04_pf_PhotonEt"          );
  produces<vector<float> >          ( branchprefix_ + "isoR04pfNeutralHadronEtHighThreshold").setBranchAlias( aliasprefix_ + "_isoR04_pf_NeutralHadronEtHighThreshold");
  produces<vector<float> >          ( branchprefix_ + "isoR04pfPhotonEtHighThreshold"       ).setBranchAlias( aliasprefix_ + "_isoR04_pf_PhotonEtHighThreshold"       );
  produces<vector<float> >          ( branchprefix_ + "isoR04pfPUPt"              ).setBranchAlias( aliasprefix_ + "_isoR04_pf_PUPt"              );

  produces<vector<float> >          ( branchprefix_ + "iso03pf"                   ).setBranchAlias( aliasprefix_ + "_iso03_pf"            ); // pf isolation in cone of 0.3
  produces<vector<float> >          ( branchprefix_ + "iso04pf"                   ).setBranchAlias( aliasprefix_ + "_iso04_pf"            ); // pf isolation in cone of 0.4

  ///////////
  // IP 3D //
  ///////////
  
  produces<vector<float> >          ( branchprefix_ + "ip3d"                      ).setBranchAlias( aliasprefix_ + "_ip3d"                ); // Ip3d from standard vertex
  produces<vector<float> >          ( branchprefix_ + "ip3derr"                   ).setBranchAlias( aliasprefix_ + "_ip3derr"             ); // Ip3d error from standard vertex

} // end Constructor

void MuonMaker::beginJob () {}  // method called once each job just before starting event loop
void MuonMaker::endJob   () {}  // method called once each job just after ending the event loop


//////////////
// Producer //
//////////////

void MuonMaker::produce(Event& iEvent, const EventSetup& iSetup) {

  /////////
  // STA //
  /////////

  auto_ptr<vector<float> >         vector_mus_sta_d0                      ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_z0                      ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_d0Err                   ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_z0Err                   ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_d0corr                  ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_z0corr                  ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_qoverp                  ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_qoverpError             ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_chi2                    ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_sta_ndof                    ( new vector<float>   );
  auto_ptr<vector<int> >           vector_mus_sta_validHits               ( new vector<int>     );
  auto_ptr<vector<LorentzVector> > vector_mus_sta_p4                      ( new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> > vector_mus_sta_vertex_p4               ( new vector<LorentzVector> );


  ////////////  
  // Global //
  ////////////

  auto_ptr<vector<float> >         vector_mus_gfit_d0                     ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_z0                     ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_d0Err                  ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_z0Err                  ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_d0corr                 ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_z0corr                 ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_qoverp                 ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_qoverpError            ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_chi2                   ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_gfit_ndof                   ( new vector<float>         );
  auto_ptr<vector<int> >           vector_mus_gfit_validHits              ( new vector<int>           );
  auto_ptr<vector<int> >           vector_mus_gfit_validSTAHits           ( new vector<int>           );
  auto_ptr<vector<int> >           vector_mus_gfit_validSiHits            ( new vector<int>           );
  auto_ptr<vector<LorentzVector> > vector_mus_gfit_p4                     ( new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_gfit_vertex_p4              ( new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_gfit_outerPos_p4            ( new vector<LorentzVector> );


  /////////////
  // Quality //
  /////////////

  auto_ptr<vector<bool> >   vector_mus_updatedSta          ( new vector<bool>   );
  auto_ptr<vector<bool> >   vector_mus_tightMatch          ( new vector<bool>   );
  auto_ptr<vector<float> >  vector_mus_trkKink             ( new vector<float>  );
  auto_ptr<vector<float> >  vector_mus_glbKink             ( new vector<float>  );
  auto_ptr<vector<float> >  vector_mus_trkRelChi2          ( new vector<float>  );
  auto_ptr<vector<float> >  vector_mus_staRelChi2          ( new vector<float>  );
  auto_ptr<vector<float> >  vector_mus_chi2LocalPosition   ( new vector<float>  );
  auto_ptr<vector<float> >  vector_mus_chi2LocalMomentum   ( new vector<float>  );
  auto_ptr<vector<float> >  vector_mus_localDistance       ( new vector<float>  );
  auto_ptr<vector<float> >  vector_mus_globalDeltaEtaPhi   ( new vector<float>  );
  auto_ptr<vector<float> >  vector_mus_glbTrackProbability ( new vector<float>  );


  ////////////
  // Timing //
  ////////////
  
  auto_ptr<vector<int> >           vector_mus_timeNumStationsUsed         ( new vector<int>     );
  auto_ptr<vector<float> >         vector_mus_timeAtIpInOut               ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_timeAtIpInOutErr            ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_timeAtIpOutIn               ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_timeAtIpOutInErr            ( new vector<float>   );
  auto_ptr<vector<int> >           vector_mus_timeDirection               ( new vector<int>     );


  ///////////////////////////
  // Cosmics Compatibility //
  ///////////////////////////

  auto_ptr<vector<float> > cosmicCompat     ( new vector<float> );
  auto_ptr<vector<float> > timeCompat       ( new vector<float> );
  auto_ptr<vector<float> > backToBackCompat ( new vector<float> );
  auto_ptr<vector<float> > overlapCompat    ( new vector<float> );
  auto_ptr<vector<float> > vertexCompat     ( new vector<float> );


  ///////////
  // Muons //
  ///////////

  auto_ptr<vector<int> >           vector_mus_type                    ( new vector<int>           );        
  auto_ptr<vector<int> >           vector_mus_goodmask                ( new vector<int>           );        
  auto_ptr<vector<int> >           vector_mus_nmatches                ( new vector<int>           );
  auto_ptr<vector<int> >           vector_mus_charge                  ( new vector<int>           );        
  auto_ptr<vector<int> >           vector_mus_nOverlaps               ( new vector<int>           );
  auto_ptr<vector<int> >           vector_mus_overlap0                ( new vector<int>           );
  auto_ptr<vector<int> >           vector_mus_overlap1                ( new vector<int>           );
  auto_ptr<vector<float> >         vector_mus_caloCompatibility       ( new vector<float>         );
  auto_ptr<vector<float> >         vector_mus_segmCompatibility       ( new vector<float>         );
  auto_ptr<vector<LorentzVector> > vector_mus_p4                      ( new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_ecalpos_p4              ( new vector<LorentzVector> );
  auto_ptr<vector<int> >           vector_mus_numberOfMatchedStations ( new vector<int>           );

  ////////
  // ID //
  ////////

  auto_ptr<vector<int> >           vector_mus_pid_TMLastStationLoose      ( new vector<int>     );
  auto_ptr<vector<int> >           vector_mus_pid_TMLastStationTight      ( new vector<int>     );
  auto_ptr<vector<int> >           vector_mus_pid_TM2DCompatibilityLoose  ( new vector<int>     );
  auto_ptr<vector<int> >           vector_mus_pid_TM2DCompatibilityTight  ( new vector<int>     );
  auto_ptr<vector<int> >           vector_mus_pid_PFMuon                  ( new vector<int>     );


  ////////////
  // Energy //
  ////////////

  auto_ptr<vector<float> >         vector_mus_e_em                ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_had               ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_ho                ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_emS9              ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_hadS9             ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_e_hoS9              ( new vector<float>          );


  ///////////////
  // Isolation //
  ///////////////

  auto_ptr<vector<float> >         vector_mus_iso_trckvetoDep     ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso_ecalvetoDep     ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso_hcalvetoDep     ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso_hovetoDep       ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso03_sumPt         ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso03_emEt          ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso03_hadEt         ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso03_hoEt          ( new vector<float>          );
  auto_ptr<vector<int> >           vector_mus_iso03_ntrk          ( new vector<int>            );
  auto_ptr<vector<float> >         vector_mus_iso05_sumPt         ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso05_emEt          ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso05_hadEt         ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso05_hoEt          ( new vector<float>          );
  auto_ptr<vector<int> >           vector_mus_iso05_ntrk          ( new vector<int>            );


  ////////////
  // Tracks //
  ////////////

  auto_ptr<vector<int> >           vector_mus_muonBestTrackType   ( new vector<int>            );        
  auto_ptr<vector<LorentzVector> > vector_mus_trk_p4              ( new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> > vector_mus_vertex_p4           ( new vector<LorentzVector> );
  auto_ptr<vector<int>   >         vector_mus_trkidx              ( new vector<int>            );
  auto_ptr<vector<float> >         vector_mus_d0                  ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_z0                  ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_d0corr              ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_z0corr              ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_vertexphi           ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_chi2                ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_ndof                ( new vector<float>          );      
  auto_ptr<vector<int> >           vector_mus_validHits           ( new vector<int>            );        
  auto_ptr<vector<int> >           vector_mus_lostHits            ( new vector<int>            );        
  auto_ptr<vector<float> >         vector_mus_d0Err               ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_z0Err               ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_ptErr               ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_etaErr              ( new vector<float>          );      
  auto_ptr<vector<float> >         vector_mus_phiErr              ( new vector<float>          );      
  auto_ptr<vector<int> >           vector_mus_trk_charge          ( new vector<int>            );   
  auto_ptr<vector<float> >         vector_mus_qoverp              ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_qoverpError         ( new vector<float>          );

  auto_ptr<vector<LorentzVector> > vector_mus_fittpfms_p4                 ( new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_fitpicky_p4                 ( new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> > vector_mus_fittev_p4                   ( new vector<LorentzVector> );


  ////////
  // PF //
  ////////

  auto_ptr< vector< int> >           vector_mus_pfcharge              ( new vector<int>   );
  auto_ptr< vector< int> >           vector_mus_pfparticleId          ( new vector<int>   );
  auto_ptr< vector< int> >           vector_mus_pfflag                ( new vector<int>   );
  auto_ptr< vector< float> >         vector_mus_pfecalE               ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfhcalE               ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfrawEcalE            ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfrawHcalE            ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfpS1E                ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfpS2E                ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfdeltaP              ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfmvaepi              ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfmvaemu              ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfmvapimu             ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfmvanothinggamma     ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfmvanothingnh        ( new vector<float> );
  auto_ptr< vector< LorentzVector> > vector_mus_pfp4                  ( new vector<LorentzVector> );
  auto_ptr< vector< LorentzVector> > vector_mus_pfposAtEcalp4         ( new vector<LorentzVector> );
  /*
  auto_ptr< vector< float> >         vector_mus_pfiso03ChargedHadrons ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfiso03NeutralHadrons ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfiso03Photons        ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfiso04ChargedHadrons ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfiso04NeutralHadrons ( new vector<float> );
  auto_ptr< vector< float> >         vector_mus_pfiso04Photons        ( new vector<float> );
  */
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_ChargedHadronPt                 ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_ChargedParticlePt               ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_NeutralHadronEt                 ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_PhotonEt                        ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_sumNeutralHadronEtHighThreshold ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_sumPhotonEtHighThreshold        ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR03_pf_PUPt                            ( new vector<float>   );


  auto_ptr< vector<float> >         vector_mus_isoR04_pf_ChargedHadronPt                 ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR04_pf_ChargedParticlePt               ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR04_pf_NeutralHadronEt                 ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR04_pf_PhotonEt                        ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR04_pf_sumNeutralHadronEtHighThreshold ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR04_pf_sumPhotonEtHighThreshold        ( new vector<float>   );
  auto_ptr< vector<float> >         vector_mus_isoR04_pf_PUPt                            ( new vector<float>   );

  auto_ptr<vector<float> >         vector_mus_iso03_pf            ( new vector<float>          );
  auto_ptr<vector<float> >         vector_mus_iso04_pf            ( new vector<float>          );

  ///////////
  // IP 3D //
  ///////////

  auto_ptr<vector<float> >         vector_mus_ip3d                        ( new vector<float>   );
  auto_ptr<vector<float> >         vector_mus_ip3derr                     ( new vector<float>   );




////////////////////////////
// --- Fill Muon Data --- //
////////////////////////////

  ///////////////
  // Get Muons //
  ///////////////

  Handle<View<Muon> > muon_h;
  iEvent.getByLabel( muonsInputTag , muon_h );

  
  ////////////////////////////////
  // Get pf-Muon from reco-muon //
  ////////////////////////////////

  Handle<ValueMap<PFCandidatePtr> > pfMap;
  iEvent.getByLabel( "particleFlow"  , muonsInputTag.label(), pfMap );


  //////////////////
  // Get Vertices //
  //////////////////

  iEvent.getByLabel( vtxInputTag , vertexHandle );  

  ///////////////////////////
  // TransientTrackBuilder //
  ///////////////////////////
  ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

  ///////////////////////
  // Get PF Candidates //
  ///////////////////////

  iEvent.getByLabel( pfCandsInputTag , pfCand_h );


  /////////////////////////////////////
  // Get BeamSpot from BeamSpotMaker //
  /////////////////////////////////////

  InputTag beamSpot_tag( beamSpotInputTag.label(), "evtbsp4" );
  Handle<LorentzVector> beamSpotH;
  iEvent.getByLabel( beamSpot_tag, beamSpotH );
  const Point beamSpot = beamSpotH.isValid() ? Point(beamSpotH->x(), beamSpotH->y(), beamSpotH->z()) : Point(0,0,0);

  
  ////////////////////////// 
  // Cosmic Compatibility //
  //////////////////////////

  Handle<MuonCollection> muons;
  iEvent.getByLabel( "muons", muons );
  Handle<ValueMap<MuonCosmicCompatibility> > CosmicMap;
  iEvent.getByLabel( src_, CosmicMap );


  ///////////
  // Muons // 
  ///////////
  
  unsigned int muonIndex = 0;
  View<Muon>::const_iterator muons_end = muon_h->end();  // Iterator
  for ( View<Muon>::const_iterator muon = muon_h->begin(); muon != muons_end; ++muon ) {

    // References
    const RefToBase<Muon>         muonRef                 = muon_h->refAt(muonIndex); 
    const TrackRef                globalTrack             = muon->globalTrack();
    const TrackRef                siTrack                 = muon->innerTrack();
    const TrackRef                staTrack                = muon->outerTrack();
    const MuonQuality             quality                 = muon->combinedQuality();
    const MuonCosmicCompatibility muonCosmicCompatibility = CosmicMap.isValid() ? (*CosmicMap)[muonRef] : MuonCosmicCompatibility(); //protection for fast sim
    const VertexCollection*       vertexCollection        = vertexHandle.product();

    // Iterators
    VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();


    /////////
    // STA //
    /////////

    vector_mus_sta_d0            -> push_back( staTrack.isNonnull()  ? staTrack->d0()                   :  -9999.        );
    vector_mus_sta_z0            -> push_back( staTrack.isNonnull()  ? staTrack->dz()                   :  -9999.        );
    vector_mus_sta_d0Err         -> push_back( staTrack.isNonnull()  ? staTrack->d0Error()              :  -9999.        );
    vector_mus_sta_z0Err         -> push_back( staTrack.isNonnull()  ? staTrack->dzError()              :  -9999.        );
    vector_mus_sta_d0corr        -> push_back( staTrack.isNonnull()  ? -1*(staTrack->dxy(beamSpot))     :  -9999.        );
    vector_mus_sta_z0corr        -> push_back( staTrack.isNonnull()  ? staTrack->dz(beamSpot)           :  -9999.        );
    vector_mus_sta_qoverp        -> push_back( staTrack.isNonnull()  ? staTrack->qoverp()               :  -9999.        );
    vector_mus_sta_qoverpError   -> push_back( staTrack.isNonnull()  ? staTrack->qoverpError()          :  -9999.        );
    vector_mus_sta_chi2          -> push_back( staTrack.isNonnull()  ? staTrack->chi2()                 :  -9999.        );
    vector_mus_sta_ndof          -> push_back( staTrack.isNonnull()  ? staTrack->ndof()                 :  -9999         );
    vector_mus_sta_validHits     -> push_back( staTrack.isNonnull()  ? staTrack->numberOfValidHits()    :  -9999         );
    vector_mus_sta_p4            -> push_back( staTrack.isNonnull()  ? LorentzVector( staTrack->px() , staTrack->py() , staTrack->pz() , staTrack->p() ) : LorentzVector(0, 0, 0, 0) );
    vector_mus_sta_vertex_p4     -> push_back( staTrack.isNonnull()  ? LorentzVector( staTrack->vx() , staTrack->vy() , staTrack->vz() ,           0.0 ) : LorentzVector(0, 0, 0, 0) );


    ////////////
    // Global //
    ////////////

    vector_mus_gfit_d0           -> push_back( globalTrack.isNonnull() ? globalTrack->d0()                 :  -9999.        );
    vector_mus_gfit_z0           -> push_back( globalTrack.isNonnull() ? globalTrack->dz()                 :  -9999.        );
    vector_mus_gfit_d0Err        -> push_back( globalTrack.isNonnull() ? globalTrack->d0Error()            :  -9999.        );
    vector_mus_gfit_z0Err        -> push_back( globalTrack.isNonnull() ? globalTrack->dzError()            :  -9999.        );
    vector_mus_gfit_d0corr       -> push_back( globalTrack.isNonnull() ? -1*(globalTrack->dxy(beamSpot))   :  -9999.        );
    vector_mus_gfit_z0corr       -> push_back( globalTrack.isNonnull() ? globalTrack->dz(beamSpot)         :  -9999.        );
    vector_mus_gfit_qoverp       -> push_back( globalTrack.isNonnull() ? globalTrack->qoverp()             :  -9999.        );
    vector_mus_gfit_qoverpError  -> push_back( globalTrack.isNonnull() ? globalTrack->qoverpError()        :  -9999.        );
    vector_mus_gfit_chi2         -> push_back( globalTrack.isNonnull() ? globalTrack->chi2()               :  -9999.        );
    vector_mus_gfit_ndof         -> push_back( globalTrack.isNonnull() ? globalTrack->ndof()               :  -9999         );
    vector_mus_gfit_validHits    -> push_back( globalTrack.isNonnull() ? globalTrack->numberOfValidHits()  :  -9999         );
    vector_mus_gfit_validSTAHits -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfValidMuonHits()    : -9999         );
    vector_mus_gfit_validSiHits  -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfValidTrackerHits() : -9999         );
    vector_mus_gfit_p4           -> push_back( globalTrack.isNonnull() ? LorentzVector( globalTrack->px()                , globalTrack->py()                , globalTrack->pz()                , globalTrack->p() ) : LorentzVector(     0.0,     0.0,    0.0,      0.0) );
    vector_mus_gfit_vertex_p4    -> push_back( globalTrack.isNonnull() ? LorentzVector( globalTrack->vx()                , globalTrack->vy()                , globalTrack->vz()                ,              0.0 ) : LorentzVector( -9999.0, -9999.0, -9999.0, -9999.0) );
    vector_mus_gfit_outerPos_p4  -> push_back( globalTrack.isNonnull() ? LorentzVector( globalTrack->outerPosition().x() , globalTrack->outerPosition().y() , globalTrack->outerPosition().z() ,              0.0 ) : LorentzVector( -9999.0, -9999.0, -9999.0, -9999.0) );


    //////////////////
    // Muon Quality //
    //////////////////

    vector_mus_updatedSta          -> push_back( quality.updatedSta          );
    vector_mus_tightMatch          -> push_back( quality.tightMatch          );
    vector_mus_trkKink             -> push_back( quality.trkKink             );
    vector_mus_glbKink             -> push_back( quality.glbKink             );
    vector_mus_trkRelChi2          -> push_back( quality.trkRelChi2          );
    vector_mus_staRelChi2          -> push_back( quality.staRelChi2          );
    vector_mus_chi2LocalPosition   -> push_back( quality.chi2LocalPosition   );
    vector_mus_chi2LocalMomentum   -> push_back( quality.chi2LocalMomentum   );
    vector_mus_localDistance       -> push_back( quality.localDistance       );
    vector_mus_globalDeltaEtaPhi   -> push_back( quality.globalDeltaEtaPhi   );
    vector_mus_glbTrackProbability -> push_back( quality.glbTrackProbability );


    ////////////
    // Timing //
    ////////////

    bool timeIsValid = muon->isTimeValid();

    vector_mus_timeNumStationsUsed -> push_back( timeIsValid ? muon->time().nDof             : -9999. );
    vector_mus_timeAtIpInOut       -> push_back( timeIsValid ? muon->time().timeAtIpInOut    : -9999. );
    vector_mus_timeAtIpInOutErr    -> push_back( timeIsValid ? muon->time().timeAtIpInOutErr : -9999. );
    vector_mus_timeAtIpOutIn       -> push_back( timeIsValid ? muon->time().timeAtIpOutIn    : -9999. );
    vector_mus_timeAtIpOutInErr    -> push_back( timeIsValid ? muon->time().timeAtIpOutInErr : -9999. );
    vector_mus_timeDirection       -> push_back( timeIsValid ? muon->time().direction()      : -9999. );


    //////////////////////////
    // Cosmic Compatibility //
    //////////////////////////

    cosmicCompat     -> push_back( muonCosmicCompatibility.cosmicCompatibility     );
    timeCompat       -> push_back( muonCosmicCompatibility.timeCompatibility       );
    backToBackCompat -> push_back( muonCosmicCompatibility.backToBackCompatibility );
    overlapCompat    -> push_back( muonCosmicCompatibility.overlapCompatibility    );
    vertexCompat     -> push_back( muonCosmicCompatibility.vertexCompatibility     );


    //////////
    // Muon //
    //////////

    // Calculate Overlaps
    int mus_overlap0 = -1, mus_overlap1 = -1, muInd = -1, mus_nOverlaps = 0;
    for ( View<Muon>::const_iterator muonJ = muon_h->begin(); muonJ != muons_end; ++muonJ ) {
      muInd++;
      if ( muonJ != muon ){
        if ( muon::overlap( *muon, *muonJ ) ) {
          if ( mus_overlap0 == -1) mus_overlap0 = muInd;
          if ( mus_overlap0 != -1) mus_overlap1 = muInd;
          mus_nOverlaps++;
        }
      }
    }

    // Calculate Muon position at ECAL
    math::XYZPoint ecal_p( -9999.0, -9999.0, -9999.0 );
    if( muon->isEnergyValid() ) ecal_p = muon->calEnergy().ecal_position;

    // Calculate Mask
    int goodMask = 0;
    for ( int iG = 0; iG < 24; ++iG ) { //overkill here
      if( isGoodMuon( *muon, (muon::SelectionType)iG ) ) goodMask |= (1<<iG);
    }

    
    /////////////////////
    // Muon Quantities //
    /////////////////////

    vector_mus_type                    -> push_back( muon->type()                                              );
    vector_mus_goodmask                -> push_back( goodMask                                                  );
    vector_mus_charge                  -> push_back( muon->charge()                                            );
    vector_mus_nmatches                -> push_back( muon->isMatchesValid() ? muon->numberOfMatches() :  -9999 );
    vector_mus_nOverlaps               -> push_back( mus_nOverlaps                                             );
    vector_mus_overlap0                -> push_back( mus_overlap0                                              );
    vector_mus_overlap1                -> push_back( mus_overlap1                                              );
    vector_mus_caloCompatibility       -> push_back( muon->caloCompatibility()                                 );
    vector_mus_segmCompatibility       -> push_back( muon::segmentCompatibility(*muon)                         );
    vector_mus_p4                      -> push_back( LorentzVector( muon->p4()                              )  );
    vector_mus_ecalpos_p4              -> push_back( LorentzVector( ecal_p.x(), ecal_p.y(), ecal_p.z(), 0.0 )  );
    vector_mus_numberOfMatchedStations ->push_back( muon->numberOfMatchedStations()                            );


    ////////
    // ID //
    ////////

    bool matchIsValid = muon->isMatchesValid();

    vector_mus_pid_TMLastStationLoose     -> push_back( matchIsValid ? muon::isGoodMuon( *muon, muon::TMLastStationLoose     ) : -9999  );
    vector_mus_pid_TMLastStationTight     -> push_back( matchIsValid ? muon::isGoodMuon( *muon, muon::TMLastStationTight     ) : -9999  );
    vector_mus_pid_TM2DCompatibilityLoose -> push_back( matchIsValid ? muon::isGoodMuon( *muon, muon::TM2DCompatibilityLoose ) : -9999  );
    vector_mus_pid_TM2DCompatibilityTight -> push_back( matchIsValid ? muon::isGoodMuon( *muon, muon::TM2DCompatibilityTight ) : -9999  );
    vector_mus_pid_PFMuon                 -> push_back( muon->isPFMuon() );

    ////////////
    // Energy //
    ////////////

    bool energyIsValid = muon->isEnergyValid();

    vector_mus_e_em               -> push_back( energyIsValid ? muon->calEnergy().em                                 :  -9999.       );
    vector_mus_e_had              -> push_back( energyIsValid ? muon->calEnergy().had                                :  -9999.       );
    vector_mus_e_ho               -> push_back( energyIsValid ? muon->calEnergy().ho                                 :  -9999.       );
    vector_mus_e_emS9             -> push_back( energyIsValid ? muon->calEnergy().emS9                               :  -9999.       );
    vector_mus_e_hadS9            -> push_back( energyIsValid ? muon->calEnergy().hadS9                              :  -9999.       );
    vector_mus_e_hoS9             -> push_back( energyIsValid ? muon->calEnergy().hoS9                               :  -9999.       );


    ///////////////
    // Isolation //
    ///////////////

    vector_mus_iso_trckvetoDep    -> push_back( muon->isEnergyValid()    ? muon->isolationR03().trackerVetoPt  : -9999.        );
    vector_mus_iso_ecalvetoDep    -> push_back( muon->isEnergyValid()    ? muon->isolationR03().emVetoEt       : -9999.        );
    vector_mus_iso_hcalvetoDep    -> push_back( muon->isEnergyValid()    ? muon->isolationR03().hadVetoEt      : -9999.        );
    vector_mus_iso_hovetoDep      -> push_back( muon->isEnergyValid()    ? muon->isolationR03().hoVetoEt       : -9999.        );
    vector_mus_iso03_sumPt        -> push_back( muon->isIsolationValid() ? muon->isolationR03().sumPt          : -9999.        );
    vector_mus_iso03_emEt         -> push_back( muon->isIsolationValid() ? muon->isolationR03().emEt           : -9999.        );
    vector_mus_iso03_hadEt        -> push_back( muon->isIsolationValid() ? muon->isolationR03().hadEt          : -9999.        );
    vector_mus_iso03_hoEt         -> push_back( muon->isIsolationValid() ? muon->isolationR03().hoEt           : -9999.        );
    vector_mus_iso03_ntrk         -> push_back( muon->isIsolationValid() ? muon->isolationR03().nTracks        : -9999         );
    vector_mus_iso05_sumPt        -> push_back( muon->isIsolationValid() ? muon->isolationR05().sumPt          : -9999.        );
    vector_mus_iso05_emEt         -> push_back( muon->isIsolationValid() ? muon->isolationR05().emEt           : -9999.        );
    vector_mus_iso05_hadEt        -> push_back( muon->isIsolationValid() ? muon->isolationR05().hadEt          : -9999.        );
    vector_mus_iso05_hoEt         -> push_back( muon->isIsolationValid() ? muon->isolationR05().hoEt           : -9999.        );
    vector_mus_iso05_ntrk         -> push_back( muon->isIsolationValid() ? muon->isolationR05().nTracks        : -9999         );


    ////////////
    // Tracks //
    ////////////

    vector_mus_muonBestTrackType  -> push_back( muon->muonBestTrackType() );
    vector_mus_trk_p4             -> push_back( siTrack.isNonnull()     ? LorentzVector( siTrack.get()->px() , siTrack.get()->py() , siTrack.get()->pz() , siTrack.get()->p() ) : LorentzVector(     0.0,     0.0,     0.0,     0.0) );
    vector_mus_vertex_p4          -> push_back( siTrack.isNonnull()     ? LorentzVector( siTrack->vx()       , siTrack->vy()       , siTrack->vz()       , 0.0                ) : LorentzVector( -9999.0, -9999.0, -9999.0, -9999.0) );
    vector_mus_trkidx             -> push_back( siTrack.isNonnull()     ? static_cast<int>(siTrack.key())                      : -9999         );
    vector_mus_d0                 -> push_back( siTrack.isNonnull()     ? siTrack->d0()                                        : -9999.        );
    vector_mus_z0                 -> push_back( siTrack.isNonnull()     ? siTrack->dz()                                        : -9999.        );
    vector_mus_d0corr             -> push_back( siTrack.isNonnull()     ? -1*(siTrack->dxy(beamSpot))                          : -9999.        );
    vector_mus_z0corr             -> push_back( siTrack.isNonnull()     ? siTrack->dz(beamSpot)                                : -9999.        );
    vector_mus_vertexphi          -> push_back( siTrack.isNonnull()     ? atan2( siTrack->vy(), siTrack->vx() )                : -9999.        );
    vector_mus_chi2               -> push_back( siTrack.isNonnull()     ? siTrack->chi2()                                      : -9999.        );
    vector_mus_ndof               -> push_back( siTrack.isNonnull()     ? siTrack->ndof()                                      : -9999.        );
    vector_mus_validHits          -> push_back( siTrack.isNonnull()     ? siTrack->numberOfValidHits()                         : -9999         );
    vector_mus_lostHits           -> push_back( siTrack.isNonnull()     ? siTrack->numberOfLostHits()                          : -9999         );
    vector_mus_d0Err              -> push_back( siTrack.isNonnull()     ? siTrack->d0Error()                                   :  -9999.       );
    vector_mus_z0Err              -> push_back( siTrack.isNonnull()     ? siTrack->dzError()                                   :  -9999.       );
    vector_mus_ptErr              -> push_back( siTrack.isNonnull()     ? siTrack->ptError()                                   :  -9999.       );
    vector_mus_etaErr             -> push_back( siTrack.isNonnull()     ? siTrack->etaError()                                  :  -9999.       );
    vector_mus_phiErr             -> push_back( siTrack.isNonnull()     ? siTrack->phiError()                                  :  -9999.       );
    vector_mus_trk_charge         -> push_back( siTrack.isNonnull()     ? siTrack->charge()                                    :  -9999        );
    vector_mus_qoverp             -> push_back( siTrack.isNonnull()     ? siTrack->qoverp()                                    :  -9999.       );
    vector_mus_qoverpError        -> push_back( siTrack.isNonnull()     ? siTrack->qoverpError()                               :  -9999.       );

    //   
    if( !muon->isGlobalMuon() ) { // Muon is not global

      vector_mus_fittpfms_p4    -> push_back( LorentzVector( 0, 0, 0, 0 ) );
      vector_mus_fitpicky_p4    -> push_back( LorentzVector( 0, 0, 0, 0 ) );
      vector_mus_fittev_p4      -> push_back( LorentzVector( 0, 0, 0, 0 ) );

    }
    else {  // Muon is global

      ///////////////
      // First Hit //
      ///////////////

      TrackRef tmpRef = muon->muonTrackFromMap( Muon::TPFMS );
      if( tmpRef.isAvailable() ) { 
        vector_mus_fittpfms_p4->push_back( LorentzVector( tmpRef->px(), tmpRef->py(), tmpRef->pz(), tmpRef->p() ) );
      }
      else{
        vector_mus_fittpfms_p4->push_back( LorentzVector( 0, 0, 0, 0 ) );
      }

      ///////////
      // Picky //
      ///////////

      tmpRef = muon->muonTrackFromMap( Muon::Picky );
      if( tmpRef.isAvailable() ) {
        vector_mus_fitpicky_p4->push_back( LorentzVector( tmpRef->px(), tmpRef->py(), tmpRef->pz(), tmpRef->p() ) );
      }
      else {
        vector_mus_fitpicky_p4->push_back( LorentzVector( 0, 0, 0, 0 ) );
      }


      /////////
      // TeV //
      /////////  

      tmpRef = muon::tevOptimized( *muon ).first;
      if( tmpRef.isAvailable() ) { 
        vector_mus_fittev_p4 -> push_back( LorentzVector( tmpRef->px(), tmpRef->py(), tmpRef->pz(), tmpRef->p() ) );
      }
      else {
        vector_mus_fittev_p4 -> push_back( LorentzVector( 0, 0, 0, 0 ) );
      }

    } // end Muon is global




    ////////
    // PF //
    ////////

    // PF p4
    vector_mus_pfp4 -> push_back( LorentzVector( muon->pfP4() ) );

    // PF Isolation
    MuonPFIsolation pfStructR03 = muon->pfIsolationR03();
    MuonPFIsolation pfStructR04 = muon->pfIsolationR04();

    vector_mus_isoR03_pf_ChargedHadronPt                 -> push_back( pfStructR03.sumChargedHadronPt              );
    vector_mus_isoR03_pf_ChargedParticlePt               -> push_back( pfStructR03.sumChargedParticlePt            );
    vector_mus_isoR03_pf_NeutralHadronEt                 -> push_back( pfStructR03.sumNeutralHadronEt              );
    vector_mus_isoR03_pf_PhotonEt                        -> push_back( pfStructR03.sumPhotonEt                     );
    vector_mus_isoR03_pf_sumNeutralHadronEtHighThreshold -> push_back( pfStructR03.sumNeutralHadronEtHighThreshold );
    vector_mus_isoR03_pf_sumPhotonEtHighThreshold        -> push_back( pfStructR03.sumPhotonEtHighThreshold        );
    vector_mus_isoR03_pf_PUPt                            -> push_back( pfStructR03.sumPUPt                         );

    vector_mus_isoR04_pf_ChargedHadronPt                 -> push_back( pfStructR04.sumChargedHadronPt              );
    vector_mus_isoR04_pf_ChargedParticlePt               -> push_back( pfStructR04.sumChargedParticlePt            );
    vector_mus_isoR04_pf_NeutralHadronEt                 -> push_back( pfStructR04.sumNeutralHadronEt              );
    vector_mus_isoR04_pf_PhotonEt                        -> push_back( pfStructR04.sumPhotonEt                     );
    vector_mus_isoR04_pf_sumNeutralHadronEtHighThreshold -> push_back( pfStructR04.sumNeutralHadronEtHighThreshold );
    vector_mus_isoR04_pf_sumPhotonEtHighThreshold        -> push_back( pfStructR04.sumPhotonEtHighThreshold        );
    vector_mus_isoR04_pf_PUPt                            -> push_back( pfStructR04.sumPUPt                         );

    // Other PF
    PFCandidatePtr pfCandRef = (*pfMap)[muonRef];
    if ( pfCandRef.isNonnull() ) {

      // Flags
      int pfflags = 0;
      for( unsigned int i = 0; i < 17; i++ ) {
        if( pfCandRef->flag( (PFCandidate::Flags) i ) ) pfflags |= (1<<i);
      }

      //
      vector_mus_pfcharge              ->push_back( pfCandRef->charge()                                                     );
      vector_mus_pfparticleId          ->push_back( pfCandRef->translateTypeToPdgId( pfCandRef->particleId() )              );
      vector_mus_pfflag                ->push_back( pfflags                                                                 );
      vector_mus_pfecalE               ->push_back( isfinite( pfCandRef->ecalEnergy() ) ? pfCandRef->ecalEnergy() : -9999.0 );
      vector_mus_pfrawEcalE            ->push_back( pfCandRef->rawEcalEnergy()                                              );
      vector_mus_pfrawHcalE            ->push_back( pfCandRef->rawHcalEnergy()                                              );
      vector_mus_pfpS1E                ->push_back( pfCandRef->pS1Energy()                                                  );
      vector_mus_pfpS2E                ->push_back( pfCandRef->pS2Energy()                                                  );
      vector_mus_pfdeltaP              ->push_back( pfCandRef->deltaP()                                                     );
      vector_mus_pfmvaepi              ->push_back( pfCandRef->mva_e_pi()                                                   );
      vector_mus_pfmvaemu              ->push_back( pfCandRef->mva_e_mu()                                                   );
      vector_mus_pfmvapimu             ->push_back( pfCandRef->mva_pi_mu()                                                  );
      vector_mus_pfmvanothinggamma     ->push_back( pfCandRef->mva_nothing_gamma()                                          );
      vector_mus_pfmvanothingnh        ->push_back( pfCandRef->mva_nothing_nh()                                             );
      vector_mus_pfposAtEcalp4         ->push_back( LorentzVector( pfCandRef->positionAtECALEntrance().x() ,
                                                                   pfCandRef->positionAtECALEntrance().y() ,
                                                                   pfCandRef->positionAtECALEntrance().z() ,
                                                                   0.0 )                                                    );
      //vector_mus_pfiso03ChargedHadrons ->push_back(                                                                     ); 
      //vector_mus_pfiso03NeutralHadrons ->push_back(                                                                     ); 
      //vector_mus_pfiso03Photons        ->push_back(                                                                     ); 
      //vector_mus_pfiso04ChargedHadrons ->push_back(                                                                     ); 
      //vector_mus_pfiso04NeutralHadrons ->push_back(                                                                     ); 
      //vector_mus_pfiso04Photons        ->push_back(                                                                     ); 

    } else {

      //
      vector_mus_pfcharge              ->push_back( -9999.0 );
      vector_mus_pfparticleId          ->push_back( -9999.0 );
      vector_mus_pfflag                ->push_back( -9999.0 );
      vector_mus_pfecalE               ->push_back( -9999.0 );
      vector_mus_pfrawEcalE            ->push_back( -9999.0 );
      vector_mus_pfrawHcalE            ->push_back( -9999.0 );
      vector_mus_pfpS1E                ->push_back( -9999.0 );
      vector_mus_pfpS2E                ->push_back( -9999.0 );
      vector_mus_pfdeltaP              ->push_back( -9999.0 );
      vector_mus_pfmvaepi              ->push_back( -9999.0 );
      vector_mus_pfmvaemu              ->push_back( -9999.0 );
      vector_mus_pfmvapimu             ->push_back( -9999.0 );
      vector_mus_pfmvanothinggamma     ->push_back( -9999.0 );
      vector_mus_pfmvanothingnh        ->push_back( -9999.0 );
      vector_mus_pfposAtEcalp4         ->push_back( LorentzVector( 0.0, 0.0, 0.0, 0.0 ) ) ;
      //vector_mus_pfiso03ChargedHadrons ->push_back(                                                                     ); 
      //vector_mus_pfiso03NeutralHadrons ->push_back(                                                                     ); 
      //vector_mus_pfiso03Photons        ->push_back(                                                                     ); 
      //vector_mus_pfiso04ChargedHadrons ->push_back(                                                                     ); 
      //vector_mus_pfiso04NeutralHadrons ->push_back(                                                                     ); 
      //vector_mus_pfiso04Photons        ->push_back(                                                                     ); 

    } //


    //////////////////    
    // PF Isolation //
    //////////////////

    for ( VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx ) {
      if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0 ) {
        firstGoodVertex = vtx;
        break;
      }
    }
    if ( firstGoodVertex!=vertexCollection->end() ) {
      vector_mus_iso03_pf->push_back( muonIsoValuePF( *muon, *firstGoodVertex, 0.3, 1.0, 0.1, 0) );
      vector_mus_iso04_pf->push_back( muonIsoValuePF( *muon, *firstGoodVertex, 0.4, 1.0, 0.1, 0) );
    } else {
      vector_mus_iso03_pf->push_back( -9999. );
      vector_mus_iso04_pf->push_back( -9999. );
    }

    ///////////
    // IP 3D //
    ///////////

    if ( siTrack.isNonnull() && firstGoodVertex != vertexCollection->end() ) {
      
      TransientTrack tt       = theTTBuilder->build( siTrack );
      Measurement1D ip3D      = IPTools::absoluteImpactParameter3D( tt, *firstGoodVertex ).second;
      vector_mus_ip3d         -> push_back( ip3D.value()                    );
      vector_mus_ip3derr      -> push_back( ip3D.error()                    );

    } else {

      vector_mus_ip3d         -> push_back( -9999. );
      vector_mus_ip3derr      -> push_back( -9999. );

    } //


    //  
    muonIndex++;

  } // end loop on muons










  /////////
  // STA //
  /////////

  iEvent.put( vector_mus_sta_d0                       , branchprefix_ + "stad0"              );
  iEvent.put( vector_mus_sta_z0                       , branchprefix_ + "staz0"              );
  iEvent.put( vector_mus_sta_d0Err                    , branchprefix_ + "stad0Err"           );
  iEvent.put( vector_mus_sta_z0Err                    , branchprefix_ + "staz0Err"           );
  iEvent.put( vector_mus_sta_d0corr                   , branchprefix_ + "stad0corr"          );
  iEvent.put( vector_mus_sta_z0corr                   , branchprefix_ + "staz0corr"          );
  iEvent.put( vector_mus_sta_qoverp                   , branchprefix_ + "staqoverp"          );
  iEvent.put( vector_mus_sta_qoverpError              , branchprefix_ + "staqoverpError"     );
  iEvent.put( vector_mus_sta_chi2                     , branchprefix_ + "stachi2"            );
  iEvent.put( vector_mus_sta_ndof                     , branchprefix_ + "standof"            );
  iEvent.put( vector_mus_sta_validHits                , branchprefix_ + "stavalidHits"       );
  iEvent.put( vector_mus_sta_p4                       , branchprefix_ + "stap4"              );
  iEvent.put( vector_mus_sta_vertex_p4                , branchprefix_ + "stavertexp4"        );
  

  ////////////                                                                       
  // Global //
  /////////////

  iEvent.put( vector_mus_gfit_d0                      , branchprefix_ + "gfitd0"             );
  iEvent.put( vector_mus_gfit_z0                      , branchprefix_ + "gfitz0"             );
  iEvent.put( vector_mus_gfit_d0Err                   , branchprefix_ + "gfitd0Err"          );
  iEvent.put( vector_mus_gfit_z0Err                   , branchprefix_ + "gfitz0Err"          );
  iEvent.put( vector_mus_gfit_d0corr                  , branchprefix_ + "gfitd0corr"         );
  iEvent.put( vector_mus_gfit_z0corr                  , branchprefix_ + "gfitz0corr"         );
  iEvent.put( vector_mus_gfit_qoverp                  , branchprefix_ + "gfitqoverp"         );
  iEvent.put( vector_mus_gfit_qoverpError             , branchprefix_ + "gfitqoverpError"    );
  iEvent.put( vector_mus_gfit_chi2                    , branchprefix_ + "gfitchi2"           );
  iEvent.put( vector_mus_gfit_ndof                    , branchprefix_ + "gfitndof"           );
  iEvent.put( vector_mus_gfit_validHits               , branchprefix_ + "gfitvalidHits"      );
  iEvent.put( vector_mus_gfit_validSTAHits            , branchprefix_ + "gfitvalidSTAHits"   );
  iEvent.put( vector_mus_gfit_validSiHits             , branchprefix_ + "gfitvalidSiHits"    );
  iEvent.put( vector_mus_gfit_p4                      , branchprefix_ + "gfitp4"             );
  iEvent.put( vector_mus_gfit_vertex_p4               , branchprefix_ + "gfitvertexp4"             );
  iEvent.put( vector_mus_gfit_outerPos_p4             , branchprefix_ + "gfitouterPosp4"           );


  //////////////////
  // Muon Quality //
  //////////////////

  iEvent.put( vector_mus_updatedSta         , branchprefix_ + "updatedSta"         );
  iEvent.put( vector_mus_tightMatch         , branchprefix_ + "tightMatch"         );
  iEvent.put( vector_mus_trkKink            , branchprefix_ + "trkKink"            );
  iEvent.put( vector_mus_glbKink            , branchprefix_ + "glbKink"            );
  iEvent.put( vector_mus_trkRelChi2         , branchprefix_ + "trkRelChi2"         );
  iEvent.put( vector_mus_staRelChi2         , branchprefix_ + "staRelChi2"         );
  iEvent.put( vector_mus_chi2LocalPosition  , branchprefix_ + "chi2LocalPosition"  );
  iEvent.put( vector_mus_chi2LocalMomentum  , branchprefix_ + "chi2LocalMomentum"  );
  iEvent.put( vector_mus_localDistance      , branchprefix_ + "localDistance"      );
  iEvent.put( vector_mus_globalDeltaEtaPhi  , branchprefix_ + "globalDeltaEtaPhi"  );
  iEvent.put( vector_mus_glbTrackProbability, branchprefix_ + "glbTrackProbability");


  ////////////
  // Timing //
  ////////////

  iEvent.put( vector_mus_timeNumStationsUsed          , branchprefix_ + "timeNumStationsUsed"      );
  iEvent.put( vector_mus_timeAtIpInOut                , branchprefix_ + "timeAtIpInOut"            );
  iEvent.put( vector_mus_timeAtIpInOutErr             , branchprefix_ + "timeAtIpInOutErr"         );
  iEvent.put( vector_mus_timeAtIpOutIn                , branchprefix_ + "timeAtIpOutIn"            );
  iEvent.put( vector_mus_timeAtIpOutInErr             , branchprefix_ + "timeAtIpOutInErr"         );
  iEvent.put( vector_mus_timeDirection                , branchprefix_ + "timeDirection"            );


  //////////////////////////
  // Cosmic Compatibility //
  //////////////////////////

  iEvent.put( cosmicCompat                  , branchprefix_ + "cosmicCompat"       );
  iEvent.put( timeCompat                    , branchprefix_ + "timeCompat"         );
  iEvent.put( backToBackCompat              , branchprefix_ + "backToBackCompat"   );
  iEvent.put( overlapCompat                 , branchprefix_ + "overlapCompat"      );
  iEvent.put( vertexCompat                  , branchprefix_ + "vertexCompat"       );


  ///////////
  // Muons //
  ///////////

  iEvent.put( vector_mus_type                    , branchprefix_ + "type"                    );
  iEvent.put( vector_mus_goodmask                , branchprefix_ + "goodmask"                );
  iEvent.put( vector_mus_charge                  , branchprefix_ + "charge"                  );
  iEvent.put( vector_mus_nmatches                , branchprefix_ + "nmatches"                );
  iEvent.put( vector_mus_nOverlaps               , branchprefix_ + "nOverlaps"               );
  iEvent.put( vector_mus_overlap0                , branchprefix_ + "overlap0"                );
  iEvent.put( vector_mus_overlap1                , branchprefix_ + "overlap1"                );
  iEvent.put( vector_mus_caloCompatibility       , branchprefix_ + "caloCompatibility"       );
  iEvent.put( vector_mus_segmCompatibility       , branchprefix_ + "segmCompatibility"       );
  iEvent.put( vector_mus_p4                      , branchprefix_ + "p4"                      );
  iEvent.put( vector_mus_ecalpos_p4              , branchprefix_ + "ecalposp4"               );
  iEvent.put( vector_mus_numberOfMatchedStations , branchprefix_ + "numberOfMatchedStations" );

  ////////
  // ID //
  ////////

  iEvent.put( vector_mus_pid_TMLastStationLoose       , branchprefix_ + "pidTMLastStationLoose"    );
  iEvent.put( vector_mus_pid_TMLastStationTight       , branchprefix_ + "pidTMLastStationTight"    );
  iEvent.put( vector_mus_pid_TM2DCompatibilityLoose   , branchprefix_ + "pidTM2DCompatibilityLoose");
  iEvent.put( vector_mus_pid_TM2DCompatibilityTight   , branchprefix_ + "pidTM2DCompatibilityTight");
  iEvent.put( vector_mus_pid_PFMuon                   , branchprefix_ + "pidPFMuon");


  ////////////
  // Energy //
  ////////////

  iEvent.put( vector_mus_e_em               , branchprefix_ + "eem"                );
  iEvent.put( vector_mus_e_had              , branchprefix_ + "ehad"               );
  iEvent.put( vector_mus_e_ho               , branchprefix_ + "eho"                );
  iEvent.put( vector_mus_e_emS9             , branchprefix_ + "eemS9"              );
  iEvent.put( vector_mus_e_hadS9            , branchprefix_ + "ehadS9"             );
  iEvent.put( vector_mus_e_hoS9             , branchprefix_ + "ehoS9"              );

  
  ///////////////
  // Isolation //
  ///////////////

  iEvent.put( vector_mus_iso_trckvetoDep    , branchprefix_ + "isotrckvetoDep"     );
  iEvent.put( vector_mus_iso_ecalvetoDep    , branchprefix_ + "isoecalvetoDep"     );
  iEvent.put( vector_mus_iso_hcalvetoDep    , branchprefix_ + "isohcalvetoDep"     );
  iEvent.put( vector_mus_iso_hovetoDep      , branchprefix_ + "isohovetoDep"       );
  iEvent.put( vector_mus_iso03_sumPt        , branchprefix_ + "iso03sumPt"         );
  iEvent.put( vector_mus_iso03_emEt         , branchprefix_ + "iso03emEt"          );
  iEvent.put( vector_mus_iso03_hadEt        , branchprefix_ + "iso03hadEt"         );
  iEvent.put( vector_mus_iso03_hoEt         , branchprefix_ + "iso03hoEt"          );
  iEvent.put( vector_mus_iso03_ntrk         , branchprefix_ + "iso03ntrk"          );
  iEvent.put( vector_mus_iso05_sumPt        , branchprefix_ + "iso05sumPt"         );
  iEvent.put( vector_mus_iso05_emEt         , branchprefix_ + "iso05emEt"          );
  iEvent.put( vector_mus_iso05_hadEt        , branchprefix_ + "iso05hadEt"         );
  iEvent.put( vector_mus_iso05_hoEt         , branchprefix_ + "iso05hoEt"          );
  iEvent.put( vector_mus_iso05_ntrk         , branchprefix_ + "iso05ntrk"          );


  ////////////
  // Tracks //
  ////////////

  iEvent.put( vector_mus_muonBestTrackType  , branchprefix_ + "muonBestTrackType"  );
  iEvent.put( vector_mus_trk_p4             , branchprefix_ + "trkp4"              );
  iEvent.put( vector_mus_vertex_p4          , branchprefix_ + "vertexp4"           );
  iEvent.put( vector_mus_trkidx             , branchprefix_ + "trkidx"             );
  iEvent.put( vector_mus_d0                 , branchprefix_ + "d0"                 );
  iEvent.put( vector_mus_z0                 , branchprefix_ + "z0"                 );
  iEvent.put( vector_mus_d0corr             , branchprefix_ + "d0corr"             );
  iEvent.put( vector_mus_z0corr             , branchprefix_ + "z0corr"             );
  iEvent.put( vector_mus_vertexphi          , branchprefix_ + "vertexphi"          );
  iEvent.put( vector_mus_chi2               , branchprefix_ + "chi2"               );
  iEvent.put( vector_mus_ndof               , branchprefix_ + "ndof"               );
  iEvent.put( vector_mus_validHits          , branchprefix_ + "validHits"          );
  iEvent.put( vector_mus_lostHits           , branchprefix_ + "lostHits"           );
  iEvent.put( vector_mus_d0Err              , branchprefix_ + "d0Err"              );
  iEvent.put( vector_mus_z0Err              , branchprefix_ + "z0Err"              );
  iEvent.put( vector_mus_ptErr              , branchprefix_ + "ptErr"              );
  iEvent.put( vector_mus_etaErr             , branchprefix_ + "etaErr"             );
  iEvent.put( vector_mus_phiErr             , branchprefix_ + "phiErr"             );
  iEvent.put( vector_mus_trk_charge         , branchprefix_ + "trkcharge"          );
  iEvent.put( vector_mus_qoverp             , branchprefix_ + "qoverp"             );
  iEvent.put( vector_mus_qoverpError        , branchprefix_ + "qoverpError"        );
          
  iEvent.put( vector_mus_fittpfms_p4                  , branchprefix_ + "fittpfmsp4"               );
  iEvent.put( vector_mus_fitpicky_p4                  , branchprefix_ + "fitpickyp4"               );
  iEvent.put( vector_mus_fittev_p4                    , branchprefix_ + "fittevp4"                 );


  ////////                                  
  // PF //
  ////////

  iEvent.put( vector_mus_pfcharge              , branchprefix_ + "pfcharge"             );
  iEvent.put( vector_mus_pfparticleId          , branchprefix_ + "pfparticleId"         );
  iEvent.put( vector_mus_pfflag                , branchprefix_ + "pfflag"               );
  iEvent.put( vector_mus_pfecalE               , branchprefix_ + "pfecalE"              );
  iEvent.put( vector_mus_pfhcalE               , branchprefix_ + "pfhcalE"              );
  iEvent.put( vector_mus_pfrawEcalE            , branchprefix_ + "pfrawEcalE"           );
  iEvent.put( vector_mus_pfrawHcalE            , branchprefix_ + "pfrawHcalE"           );
  iEvent.put( vector_mus_pfpS1E                , branchprefix_ + "pfpS1E"               );
  iEvent.put( vector_mus_pfpS2E                , branchprefix_ + "pfpS2E"               );
  iEvent.put( vector_mus_pfdeltaP              , branchprefix_ + "pfdeltaP"             );
  iEvent.put( vector_mus_pfmvaepi              , branchprefix_ + "pfmvaepi"             );
  iEvent.put( vector_mus_pfmvaemu              , branchprefix_ + "pfmvaemu"             );
  iEvent.put( vector_mus_pfmvapimu             , branchprefix_ + "pfmvapimu"            );
  iEvent.put( vector_mus_pfmvanothinggamma     , branchprefix_ + "pfmvanothinggamma"    );
  iEvent.put( vector_mus_pfmvanothingnh        , branchprefix_ + "pfmvanothingnh"       );
  iEvent.put( vector_mus_pfp4                  , branchprefix_ + "pfp4"                 );
  iEvent.put( vector_mus_pfposAtEcalp4         , branchprefix_ + "pfposAtEcalp4"        );

  /*
  iEvent.put( vector_mus_pfiso03ChargedHadrons , branchprefix_ + "pfiso03ChargedHadrons");
  iEvent.put( vector_mus_pfiso03NeutralHadrons , branchprefix_ + "pfiso03NeutralHadrons");
  iEvent.put( vector_mus_pfiso03Photons        , branchprefix_ + "pfiso03Photons"       );
  iEvent.put( vector_mus_pfiso04ChargedHadrons , branchprefix_ + "pfiso04ChargedHadrons");
  iEvent.put( vector_mus_pfiso04NeutralHadrons , branchprefix_ + "pfiso04NeutralHadrons");
  iEvent.put( vector_mus_pfiso04Photons        , branchprefix_ + "pfiso04Photons"       );
  */

  iEvent.put( vector_mus_isoR03_pf_ChargedHadronPt                , branchprefix_ + "isoR03pfChargedHadronPt"             );
  iEvent.put( vector_mus_isoR03_pf_ChargedParticlePt              , branchprefix_ + "isoR03pfChargedParticlePt"           );
  iEvent.put( vector_mus_isoR03_pf_NeutralHadronEt                , branchprefix_ + "isoR03pfNeutralHadronEt"             );
  iEvent.put( vector_mus_isoR03_pf_PhotonEt                       , branchprefix_ + "isoR03pfPhotonEt"                    );
  iEvent.put( vector_mus_isoR03_pf_sumNeutralHadronEtHighThreshold, branchprefix_ + "isoR03pfNeutralHadronEtHighThreshold");
  iEvent.put( vector_mus_isoR03_pf_sumPhotonEtHighThreshold       , branchprefix_ + "isoR03pfPhotonEtHighThreshold"       );
  iEvent.put( vector_mus_isoR03_pf_PUPt                           , branchprefix_ + "isoR03pfPUPt"                        );

  iEvent.put( vector_mus_isoR04_pf_ChargedHadronPt                , branchprefix_ + "isoR04pfChargedHadronPt"             );
  iEvent.put( vector_mus_isoR04_pf_ChargedParticlePt              , branchprefix_ + "isoR04pfChargedParticlePt"           );
  iEvent.put( vector_mus_isoR04_pf_NeutralHadronEt                , branchprefix_ + "isoR04pfNeutralHadronEt"             );
  iEvent.put( vector_mus_isoR04_pf_PhotonEt                       , branchprefix_ + "isoR04pfPhotonEt"                    );
  iEvent.put( vector_mus_isoR04_pf_sumNeutralHadronEtHighThreshold, branchprefix_ + "isoR04pfNeutralHadronEtHighThreshold");
  iEvent.put( vector_mus_isoR04_pf_sumPhotonEtHighThreshold       , branchprefix_ + "isoR04pfPhotonEtHighThreshold"       );
  iEvent.put( vector_mus_isoR04_pf_PUPt                           , branchprefix_ + "isoR04pfPUPt"                        );

  iEvent.put( vector_mus_iso03_pf , branchprefix_ + "iso03pf" );
  iEvent.put( vector_mus_iso04_pf , branchprefix_ + "iso04pf" );

  ///////////
  // IP 3D //
  ///////////

  iEvent.put( vector_mus_ip3d                         , branchprefix_ + "ip3d"               );
  iEvent.put( vector_mus_ip3derr                      , branchprefix_ + "ip3derr"            );



} //


double MuonMaker::muonIsoValuePF(const Muon& mu, const Vertex& vtx, float coner, float minptn, float dzcut, int filterId){

  float pfciso = 0;
  float pfniso = 0;

  const TrackRef siTrack  = mu.innerTrack();

  float mudz = siTrack.isNonnull() ? siTrack->dz(vtx.position()) : mu.standAloneMuon()->dz(vtx.position());

  for (PFCandidateCollection::const_iterator pf=pfCand_h->begin(); pf<pfCand_h->end(); ++pf){

    float dR = deltaR(pf->eta(), pf->phi(), mu.eta(), mu.phi());
    if (dR>coner) continue;

    int pfid = abs(pf->pdgId());
    if (filterId!=0 && filterId!=pfid) continue;

    float pfpt = pf->pt();
    if (pf->charge()==0) {
      //neutrals
      if (pfpt>minptn) pfniso+=pfpt;
    } else {
      //charged
      //avoid double counting of muon itself
      const TrackRef pfTrack  = pf->trackRef();
      if (siTrack.isNonnull()  && pfTrack.isNonnull() && siTrack.key()==pfTrack.key()) continue;
      //first check electrons with gsf track
      if (abs(pf->pdgId())==11 && pf->gsfTrackRef().isNonnull()) {
        if(fabs(pf->gsfTrackRef()->dz(vtx.position()) - mudz )<dzcut) {//dz cut
          pfciso+=pfpt;
        }
        continue;//and avoid double counting
      }
      //then check anything that has a ctf track
      if (pfTrack.isNonnull()) {//charged (with a ctf track)
        if(fabs( pfTrack->dz(vtx.position()) - mudz )<dzcut) {//dz cut
          pfciso+=pfpt;
        }
      } 
    }
  } 
  return pfciso+pfniso;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonMaker);
