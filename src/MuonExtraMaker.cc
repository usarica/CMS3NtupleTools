// -*- C++ -*-
//
// Package:    MuonExtraMaker
// Class:      MuonExtraMaker
// 
/**\class MuonExtraMaker MuonExtraMaker.cc CMS2/NtupleMaker/src/MuonExtraMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonExtraMaker.cc,v 1.68 2012/07/20 01:19:39 dbarge Exp $
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
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/VectorUtil.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CMS3/NtupleMaker/interface/MuonExtraMaker.h"

#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"


//////////////
// typedefs //
//////////////

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

MuonExtraMaker::MuonExtraMaker( const ParameterSet& iConfig ) {

    /////////////////////////////
    // Branch & Alias prefixes //
    /////////////////////////////

    aliasprefix_        = iConfig.getUntrackedParameter<string>("aliasPrefix");
    branchprefix_       = aliasprefix_; if( branchprefix_.find("_") != string::npos ) branchprefix_.replace( branchprefix_.find("_"), 1, "" );


    //////////////////////
    // Input Parameters //
    //////////////////////

    muonsToken    = consumes<View<pat::Muon> >(iConfig.getParameter<InputTag> ("muonsInputTag"   ));
    beamSpotToken  = consumes<LorentzVector>(iConfig.getParameter<edm::InputTag>("beamSpotInputTag"));
    pfCandsToken  = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<InputTag> ("pfCandsInputTag" ));
    vtxToken         = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxInputTag"));
    tevMuonsName     = iConfig.getParameter<string>   ("tevMuonsName"    );

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
    produces<vector<int> >            ( branchprefix_ + "standof"                   ).setBranchAlias( aliasprefix_ + "_sta_ndof"          ); // number of degree of freedom of the STA muon fit 
    produces<vector<int> >            ( branchprefix_ + "stavalidHits"              ).setBranchAlias( aliasprefix_ + "_sta_validHits"     ); // number of valid hits of the STA muon fit 
    produces<vector<LorentzVector> >  ( branchprefix_ + "stap4"                     ).setBranchAlias( aliasprefix_ + "_sta_p4"            ); // 
    produces<vector<LorentzVector> >  ( branchprefix_ + "stavertexp4"               ).setBranchAlias( aliasprefix_ + "_sta_vertex_p4"     );
    produces<vector<float> >          ( branchprefix_ + "stad0corrPhi"              ).setBranchAlias( aliasprefix_ + "_sta_d0corrPhi"     );
    produces<vector<float> >          ( branchprefix_ + "stad0phiCov"               ).setBranchAlias( aliasprefix_ + "_sta_d0phiCov"      );
    produces<vector<int> >            ( branchprefix_ + "staqualityMask"            ).setBranchAlias( aliasprefix_ + "_sta_qualityMask"   );
    produces<vector<int> >            ( branchprefix_ + "staalgo"                   ).setBranchAlias( aliasprefix_ + "_sta_algo"          );
    produces<vector<int> >            ( branchprefix_ + "stanlayers"                ).setBranchAlias( aliasprefix_ + "_sta_nlayers"       );
    produces<vector<int> >            ( branchprefix_ + "stanlayers3D"              ).setBranchAlias( aliasprefix_ + "_sta_nlayers3D"     );
    produces<vector<int> >            ( branchprefix_ + "stanlayersLost"            ).setBranchAlias( aliasprefix_ + "_sta_nlayersLost"   );
    produces<vector<int> >            ( branchprefix_ + "stavalidPixelHits"         ).setBranchAlias( aliasprefix_ + "_sta_validPixelHits");
    produces<vector<int> >            ( branchprefix_ + "stalostPixelHits"          ).setBranchAlias( aliasprefix_ + "_sta_lostPixelHits" );
    produces<vector<int> >            ( branchprefix_ + "staexpinnerlayer"          ).setBranchAlias( aliasprefix_ + "_sta_exp_innerlayer");
    produces<vector<int> >            ( branchprefix_ + "staexpouterlayer"          ).setBranchAlias( aliasprefix_ + "_sta_exp_outerlayer");
    produces<vector<int> >            ( branchprefix_ + "stalostHits"               ).setBranchAlias( aliasprefix_ + "_sta_lostHits"      );
    produces<vector<float> >          ( branchprefix_ + "staptErr"                  ).setBranchAlias( aliasprefix_ + "_sta_ptErr"         );
    produces<vector<float> >          ( branchprefix_ + "staetaErr"                 ).setBranchAlias( aliasprefix_ + "_sta_etaErr"        );
    produces<vector<float> >          ( branchprefix_ + "staphiErr"                 ).setBranchAlias( aliasprefix_ + "_sta_phiErr"        );
    produces<vector<int> >            ( branchprefix_ + "statrkcharge"              ).setBranchAlias( aliasprefix_ + "_sta_trk_charge"    );

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
    produces<vector<int> >            ( branchprefix_ + "gfitvalidHits"             ).setBranchAlias( aliasprefix_ + "_gfit_validHits"     ); // number of valid hits of the global muon fit 
    produces<vector<int> >            ( branchprefix_ + "gfitvalidSiHits"           ).setBranchAlias( aliasprefix_ + "_gfit_validSiHits"   ); // number of hits in the Si fit that made it into the gfit
    produces<vector<LorentzVector> >  ( branchprefix_ + "gfitvertexp4"              ).setBranchAlias( aliasprefix_ + "_gfit_vertex_p4"     );
    produces<vector<float> >          ( branchprefix_ + "gfitd0corrPhi"             ).setBranchAlias( aliasprefix_ + "_gfit_d0corrPhi"     );
    produces<vector<float> >          ( branchprefix_ + "gfitd0phiCov"              ).setBranchAlias( aliasprefix_ + "_gfit_d0phiCov"      );
    produces<vector<int> >            ( branchprefix_ + "gfitqualityMask"           ).setBranchAlias( aliasprefix_ + "_gfit_qualityMask"   );
    produces<vector<int> >            ( branchprefix_ + "gfitnlayers"               ).setBranchAlias( aliasprefix_ + "_gfit_nlayers"       );
    produces<vector<int> >            ( branchprefix_ + "gfitnlayers3D"             ).setBranchAlias( aliasprefix_ + "_gfit_nlayers3D"     );
    produces<vector<int> >            ( branchprefix_ + "gfitnlayersLost"           ).setBranchAlias( aliasprefix_ + "_gfit_nlayersLost"   );
    produces<vector<int> >            ( branchprefix_ + "gfitvalidPixelHits"        ).setBranchAlias( aliasprefix_ + "_gfit_validPixelHits");
    produces<vector<int> >            ( branchprefix_ + "gfitlostPixelHits"         ).setBranchAlias( aliasprefix_ + "_gfit_lostPixelHits" );
    produces<vector<int> >            ( branchprefix_ + "gfitexpinnerlayer"         ).setBranchAlias( aliasprefix_ + "_gfit_exp_innerlayer");
    produces<vector<int> >            ( branchprefix_ + "gfitexpouterlayer"         ).setBranchAlias( aliasprefix_ + "_gfit_exp_outerlayer");
    produces<vector<int> >            ( branchprefix_ + "gfitlostHits"              ).setBranchAlias( aliasprefix_ + "_gfit_lostHits"      );
    produces<vector<float> >          ( branchprefix_ + "gfitetaErr"                ).setBranchAlias( aliasprefix_ + "_gfit_etaErr"        );
    produces<vector<float> >          ( branchprefix_ + "gfitphiErr"                ).setBranchAlias( aliasprefix_ + "_gfit_phiErr"        );
    produces<vector<int> >            ( branchprefix_ + "gfittrkcharge"             ).setBranchAlias( aliasprefix_ + "_gfit_trk_charge"    );

  
    ////////////
    // Best   //
    ////////////

    produces<vector<float> >          ( branchprefix_ + "bfitd0"                    ).setBranchAlias( aliasprefix_ + "_bfit_d0"            ); // d0 from global fit, if it exists
    produces<vector<float> >          ( branchprefix_ + "bfitz0"                    ).setBranchAlias( aliasprefix_ + "_bfit_z0"            ); // z0 from global fit, if it exists
    produces<vector<float> >          ( branchprefix_ + "bfitd0Err"                 ).setBranchAlias( aliasprefix_ + "_bfit_d0Err"         ); // d0Err from global fit, if it exists
    produces<vector<float> >          ( branchprefix_ + "bfitz0Err"                 ).setBranchAlias( aliasprefix_ + "_bfit_z0Err"         ); // z0Err from global fit, if it exists
    produces<vector<float> >          ( branchprefix_ + "bfitd0corr"                ).setBranchAlias( aliasprefix_ + "_bfit_d0corr"        ); // Beamspot corrected d0 from global fit, if it exists
    produces<vector<float> >          ( branchprefix_ + "bfitz0corr"                ).setBranchAlias( aliasprefix_ + "_bfit_z0corr"        ); // Beamspot corrected z0 from global fit, if it exists
    produces<vector<float> >          ( branchprefix_ + "bfitqoverp"                ).setBranchAlias( aliasprefix_ + "_bfit_qoverp"        ); // global track qoverp
    produces<vector<float> >          ( branchprefix_ + "bfitqoverpError"           ).setBranchAlias( aliasprefix_ + "_bfit_qoverpError"   ); // global track qoverp error  
    produces<vector<float> >          ( branchprefix_ + "bfitchi2"                  ).setBranchAlias( aliasprefix_ + "_bfit_chi2"          ); // chi2 of the global muon fit 
    produces<vector<int> >            ( branchprefix_ + "bfitndof"                  ).setBranchAlias( aliasprefix_ + "_bfit_ndof"          ); // number of degree of freedom of the global muon fit 
    produces<vector<int> >            ( branchprefix_ + "bfitvalidHits"             ).setBranchAlias( aliasprefix_ + "_bfit_validHits"     ); // number of valid hits of the global muon fit 
    produces<vector<int> >            ( branchprefix_ + "bfitvalidSTAHits"          ).setBranchAlias( aliasprefix_ + "_bfit_validSTAHits"  ); // number of hits in the stand alone fit that made it into the bfit
    produces<vector<int> >            ( branchprefix_ + "bfitvalidSiHits"           ).setBranchAlias( aliasprefix_ + "_bfit_validSiHits"   ); // number of hits in the Si fit that made it into the bfit
    produces<vector<LorentzVector> >  ( branchprefix_ + "bfitvertexp4"              ).setBranchAlias( aliasprefix_ + "_bfit_vertex_p4"     );
    produces<vector<float> >          ( branchprefix_ + "bfitd0corrPhi"             ).setBranchAlias( aliasprefix_ + "_bfit_d0corrPhi"     );
    produces<vector<float> >          ( branchprefix_ + "bfitd0phiCov"              ).setBranchAlias( aliasprefix_ + "_bfit_d0phiCov"      );
    produces<vector<int> >            ( branchprefix_ + "bfitqualityMask"           ).setBranchAlias( aliasprefix_ + "_bfit_qualityMask"   );
    produces<vector<int> >            ( branchprefix_ + "bfitnlayers"               ).setBranchAlias( aliasprefix_ + "_bfit_nlayers"       );
    produces<vector<int> >            ( branchprefix_ + "bfitnlayers3D"             ).setBranchAlias( aliasprefix_ + "_bfit_nlayers3D"     );
    produces<vector<int> >            ( branchprefix_ + "bfitnlayersLost"           ).setBranchAlias( aliasprefix_ + "_bfit_nlayersLost"   );
    produces<vector<int> >            ( branchprefix_ + "bfitvalidPixelHits"        ).setBranchAlias( aliasprefix_ + "_bfit_validPixelHits");
    produces<vector<int> >            ( branchprefix_ + "bfitlostPixelHits"         ).setBranchAlias( aliasprefix_ + "_bfit_lostPixelHits" );
    produces<vector<int> >            ( branchprefix_ + "bfitexpinnerlayer"         ).setBranchAlias( aliasprefix_ + "_bfit_exp_innerlayer");
    produces<vector<int> >            ( branchprefix_ + "bfitexpouterlayer"         ).setBranchAlias( aliasprefix_ + "_bfit_exp_outerlayer");
    produces<vector<int> >            ( branchprefix_ + "bfitlostHits"              ).setBranchAlias( aliasprefix_ + "_bfit_lostHits"      );
    produces<vector<float> >          ( branchprefix_ + "bfitetaErr"                ).setBranchAlias( aliasprefix_ + "_bfit_etaErr"        );
    produces<vector<float> >          ( branchprefix_ + "bfitphiErr"                ).setBranchAlias( aliasprefix_ + "_bfit_phiErr"        );
    produces<vector<int> >            ( branchprefix_ + "bfittrkcharge"             ).setBranchAlias( aliasprefix_ + "_bfit_trk_charge"    );

    /////////////
    // Quality //
    /////////////

    produces<vector<bool> >           ( branchprefix_ + "updatedSta"                ).setBranchAlias( aliasprefix_ + "_updatedSta"          );  // Muon Quality - updatedSta
    produces<vector<bool> >           ( branchprefix_ + "tightMatch"                ).setBranchAlias( aliasprefix_ + "_tightMatch"          );  // Muon Quality - tightMatch
    produces<vector<float> >          ( branchprefix_ + "glbKink"                   ).setBranchAlias( aliasprefix_ + "_glbKink"             );  // Muon Quality - glbKink
    produces<vector<float> >          ( branchprefix_ + "trkRelChi2"                ).setBranchAlias( aliasprefix_ + "_trkRelChi2"          );  // Muon Quality - trkRelChi2
    produces<vector<float> >          ( branchprefix_ + "staRelChi2"                ).setBranchAlias( aliasprefix_ + "_staRelChi2"          );  // Muon Quality - staRelChi2
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

    //////////
    // Muon //
    //////////

    produces<vector<int> >            ( branchprefix_ + "goodmask"                  ).setBranchAlias( aliasprefix_ + "_goodmask"               ); // good mask
    produces<vector<int> >            ( branchprefix_ + "nmatches"                  ).setBranchAlias( aliasprefix_ + "_nmatches"               ); // number of stations with matched segments 
    produces<vector<int> >            ( branchprefix_ + "nOverlaps"                 ).setBranchAlias( aliasprefix_ + "_nOverlaps"              ); //overlap index (-1 if none)
    produces<vector<int> >            ( branchprefix_ + "overlap0"                  ).setBranchAlias( aliasprefix_ + "_overlap0"               );
    produces<vector<int> >            ( branchprefix_ + "overlap1"                  ).setBranchAlias( aliasprefix_ + "_overlap1"               );
    produces<vector<float> >          ( branchprefix_ + "mass"                      ).setBranchAlias( aliasprefix_ + "_mass"                   ); 
    produces<vector<bool> >           ( branchprefix_ + "isRPCMuon"                 ).setBranchAlias( aliasprefix_ + "_isRPCMuon"              ); 


    ////////////
    // Energy //
    ////////////

    produces<vector<float> >          ( branchprefix_ + "eem"                       ).setBranchAlias( aliasprefix_ + "_e_em"                ); // energy in crossed ECAL crystalls 
    produces<vector<float> >          ( branchprefix_ + "ehad"                      ).setBranchAlias( aliasprefix_ + "_e_had"               ); // energy in crossed HCAL towers 
    produces<vector<float> >          ( branchprefix_ + "eho"                       ).setBranchAlias( aliasprefix_ + "_e_ho"                ); // energy in crossed HO towers 
    produces<vector<float> >          ( branchprefix_ + "eemS9"                     ).setBranchAlias( aliasprefix_ + "_e_emS9"              ); // energy in 3x3 ECAL crystall shape 
    produces<vector<float> >          ( branchprefix_ + "ehadS9"                    ).setBranchAlias( aliasprefix_ + "_e_hadS9"             ); // energy in 3x3 HCAL towers 
    produces<vector<float> >          ( branchprefix_ + "ehoS9"                     ).setBranchAlias( aliasprefix_ + "_e_hoS9"              ); // energy in 3x3 HO towers 
    produces<vector<float> >          ( branchprefix_ + "emS25"         ).setBranchAlias( aliasprefix_ + "_emS25"           ); 
    produces<vector<float> >          ( branchprefix_ + "emMax"         ).setBranchAlias( aliasprefix_ + "_emMax"           ); 
    produces<vector<float> >          ( branchprefix_ + "hadMax"        ).setBranchAlias( aliasprefix_ + "_hadMax"          ); 
    produces<vector<int> >            ( branchprefix_ + "ecalrawId"     ).setBranchAlias( aliasprefix_ + "_ecal_rawId"      ); 
    produces<vector<int> >            ( branchprefix_ + "hcalrawId"     ).setBranchAlias( aliasprefix_ + "_hcal_rawId"      ); 


    ///////////////
    // Isolation //
    ///////////////
                                    
    produces<vector<float> >          ( branchprefix_ + "iso03hoEt"                 ).setBranchAlias( aliasprefix_ + "_iso03_hoEt"          ); // sum of ho Et for cone of 0.3 
    produces<vector<int> >            ( branchprefix_ + "iso03njets"                ).setBranchAlias( aliasprefix_ + "_iso03_njets"         ); // number of jets in the cone of 0.3 
    produces<vector<float> >          ( branchprefix_ + "iso05sumPt"                ).setBranchAlias( aliasprefix_ + "_iso05_sumPt"         ); // sum of track Pt for cone of 0.5 
    produces<vector<float> >          ( branchprefix_ + "iso05emEt"                 ).setBranchAlias( aliasprefix_ + "_iso05_emEt"          ); // sum of ecal Et for cone of 0.5 
    produces<vector<float> >          ( branchprefix_ + "iso05hadEt"                ).setBranchAlias( aliasprefix_ + "_iso05_hadEt"         ); // sum of hcal Et for cone of 0.5 
    produces<vector<float> >          ( branchprefix_ + "iso05hoEt"                 ).setBranchAlias( aliasprefix_ + "_iso05_hoEt"          ); // sum of ho Et for cone of 0.5 
    produces<vector<int> >            ( branchprefix_ + "iso05ntrk"                 ).setBranchAlias( aliasprefix_ + "_iso05_ntrk"          ); // number of tracks in the cone of 0.5 

    ////////////
    // Tracks //
    ////////////

    produces<vector<int> >           ( branchprefix_ + "muonBestTrackType" ).setBranchAlias( aliasprefix_ + "_muonBestTrackType" ); 
    produces<vector<LorentzVector> > ( branchprefix_ + "vertexp4"          ).setBranchAlias( aliasprefix_ + "_vertex_p4"         ); // from the silicon fit
    produces<vector<float> >         ( branchprefix_ + "d0"                ).setBranchAlias( aliasprefix_ + "_d0"                ); // impact parameter at the point of closest approach  using the tracker fit
    produces<vector<float> >         ( branchprefix_ + "z0"                ).setBranchAlias( aliasprefix_ + "_z0"                ); // z position of the point of closest approach. From the si track    
    produces<vector<float> >         ( branchprefix_ + "d0corr"            ).setBranchAlias( aliasprefix_ + "_d0corr"            ); // corrected impact parameter at the point of closest approach. From si track  
    produces<vector<float> >         ( branchprefix_ + "z0corr"            ).setBranchAlias( aliasprefix_ + "_z0corr"            ); // corrected z position of the point of closest approach. From si track    
    produces<vector<float> >         ( branchprefix_ + "vertexphi"         ).setBranchAlias( aliasprefix_ + "_vertexphi"         ); // phi angle of the point of closest approach. From si track    
    produces<vector<float> >         ( branchprefix_ + "chi2"              ).setBranchAlias( aliasprefix_ + "_chi2"              ); // chi2 of the silicon tracker fit      
    produces<vector<float> >         ( branchprefix_ + "ndof"              ).setBranchAlias( aliasprefix_ + "_ndof"              ); // number of degrees of freedom of the si tracker fit    
    produces<vector<float> >         ( branchprefix_ + "etaErr"            ).setBranchAlias( aliasprefix_ + "_etaErr"            ); // si track eta error          
    produces<vector<float> >         ( branchprefix_ + "phiErr"            ).setBranchAlias( aliasprefix_ + "_phiErr"            ); // si track phi error          
    produces<vector<int> >           ( branchprefix_ + "trkcharge"         ).setBranchAlias( aliasprefix_ + "_trk_charge"        ); // si track charge
    produces<vector<float> >         ( branchprefix_ + "qoverp"            ).setBranchAlias( aliasprefix_ + "_qoverp"            ); // si track qoverp
    produces<vector<float> >         ( branchprefix_ + "qoverpError"       ).setBranchAlias( aliasprefix_ + "_qoverpError"       ); // si track qoverp error
    produces<vector<float> >         ( branchprefix_ + "d0corrPhi"         ).setBranchAlias( aliasprefix_ + "_d0corrPhi"         );
    produces<vector<float> >         ( branchprefix_ + "d0phiCov"          ).setBranchAlias( aliasprefix_ + "_d0phiCov"          );
    produces<vector<int> >           ( branchprefix_ + "qualityMask"       ).setBranchAlias( aliasprefix_ + "_qualityMask"       );
    produces<vector<int> >           ( branchprefix_ + "nlayers3D"         ).setBranchAlias( aliasprefix_ + "_nlayers3D"         );
    produces<vector<int> >           ( branchprefix_ + "nlayersLost"       ).setBranchAlias( aliasprefix_ + "_nlayersLost"       );
    produces<vector<int> >           ( branchprefix_ + "lostPixelHits"     ).setBranchAlias( aliasprefix_ + "_lostPixelHits"     );

  
    ////////
    // PF //
    ////////

    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR03pfChargedHadronPt"   ).setBranchAlias( aliasprefix_ + "_isoMeanDRR03_pf_ChargedHadronPt"   );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR03pfChargedParticlePt" ).setBranchAlias( aliasprefix_ + "_isoMeanDRR03_pf_ChargedParticlePt" );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR03pfNeutralHadronEt"   ).setBranchAlias( aliasprefix_ + "_isoMeanDRR03_pf_NeutralHadronEt"   );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR03pfPhotonEt"          ).setBranchAlias( aliasprefix_ + "_isoMeanDRR03_pf_PhotonEt"          );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR03pfNeutralHadronEtHighThreshold").setBranchAlias( aliasprefix_ + "_isoMeanDRR03_pf_NeutralHadronEtHighThreshold");
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR03pfPhotonEtHighThreshold"       ).setBranchAlias( aliasprefix_ + "_isoMeanDRR03_pf_PhotonEtHighThreshold"       );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR03pfPUPt"              ).setBranchAlias( aliasprefix_ + "_isoMeanDRR03_pf_PUPt"              );

    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR04pfChargedHadronPt"   ).setBranchAlias( aliasprefix_ + "_isoMeanDRR04_pf_ChargedHadronPt"   );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR04pfChargedParticlePt" ).setBranchAlias( aliasprefix_ + "_isoMeanDRR04_pf_ChargedParticlePt" );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR04pfNeutralHadronEt"   ).setBranchAlias( aliasprefix_ + "_isoMeanDRR04_pf_NeutralHadronEt"   );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR04pfPhotonEt"          ).setBranchAlias( aliasprefix_ + "_isoMeanDRR04_pf_PhotonEt"          );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR04pfNeutralHadronEtHighThreshold").setBranchAlias( aliasprefix_ + "_isoMeanDRR04_pf_NeutralHadronEtHighThreshold");
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR04pfPhotonEtHighThreshold"       ).setBranchAlias( aliasprefix_ + "_isoMeanDRR04_pf_PhotonEtHighThreshold"       );
    produces<vector<float> >          ( branchprefix_ + "isoMeanDRR04pfPUPt"              ).setBranchAlias( aliasprefix_ + "_isoMeanDRR04_pf_PUPt"              );

    produces<vector<float> >          ( branchprefix_ + "isoSumDRR03pfChargedHadronPt"   ).setBranchAlias( aliasprefix_ + "_isoSumDRR03_pf_ChargedHadronPt"   );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR03pfChargedParticlePt" ).setBranchAlias( aliasprefix_ + "_isoSumDRR03_pf_ChargedParticlePt" );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR03pfNeutralHadronEt"   ).setBranchAlias( aliasprefix_ + "_isoSumDRR03_pf_NeutralHadronEt"   );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR03pfPhotonEt"          ).setBranchAlias( aliasprefix_ + "_isoSumDRR03_pf_PhotonEt"          );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR03pfNeutralHadronEtHighThreshold").setBranchAlias( aliasprefix_ + "_isoSumDRR03_pf_NeutralHadronEtHighThreshold");
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR03pfPhotonEtHighThreshold"       ).setBranchAlias( aliasprefix_ + "_isoSumDRR03_pf_PhotonEtHighThreshold"       );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR03pfPUPt"              ).setBranchAlias( aliasprefix_ + "_isoSumDRR03_pf_PUPt"              );

    produces<vector<float> >          ( branchprefix_ + "isoSumDRR04pfChargedHadronPt"   ).setBranchAlias( aliasprefix_ + "_isoSumDRR04_pf_ChargedHadronPt"   );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR04pfChargedParticlePt" ).setBranchAlias( aliasprefix_ + "_isoSumDRR04_pf_ChargedParticlePt" );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR04pfNeutralHadronEt"   ).setBranchAlias( aliasprefix_ + "_isoSumDRR04_pf_NeutralHadronEt"   );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR04pfPhotonEt"          ).setBranchAlias( aliasprefix_ + "_isoSumDRR04_pf_PhotonEt"          );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR04pfNeutralHadronEtHighThreshold").setBranchAlias( aliasprefix_ + "_isoSumDRR04_pf_NeutralHadronEtHighThreshold");
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR04pfPhotonEtHighThreshold"       ).setBranchAlias( aliasprefix_ + "_isoSumDRR04_pf_PhotonEtHighThreshold"       );
    produces<vector<float> >          ( branchprefix_ + "isoSumDRR04pfPUPt"              ).setBranchAlias( aliasprefix_ + "_isoSumDRR04_pf_PUPt"              );

} // end Constructor

void MuonExtraMaker::beginJob () {}  // method called once each job just before starting event loop
void MuonExtraMaker::endJob   () {}  // method called once each job just after ending the event loop


//////////////
// Producer //
//////////////

void MuonExtraMaker::produce(Event& iEvent, const EventSetup& iSetup) {

    /////////
    // STA //
    /////////

    unique_ptr<vector<float> >         vector_mus_sta_d0                      ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_z0                      ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_d0Err                   ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_z0Err                   ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_d0corr                  ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_z0corr                  ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_qoverp                  ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_qoverpError             ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_chi2                    ( new vector<float>         );
    unique_ptr<vector<int> >           vector_mus_sta_ndof                    ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_validHits               ( new vector<int>           );
    unique_ptr<vector<LorentzVector> > vector_mus_sta_p4                      ( new vector<LorentzVector> );
    unique_ptr<vector<LorentzVector> > vector_mus_sta_vertex_p4               ( new vector<LorentzVector> );
    unique_ptr<vector<float> >         vector_mus_sta_d0corrPhi               ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_d0phiCov                ( new vector<float>         );
    unique_ptr<vector<int> >           vector_mus_sta_qualityMask             ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_algo                    ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_nlayers                 ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_nlayers3D               ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_nlayersLost             ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_validPixelHits          ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_lostPixelHits           ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_exp_innerlayers         ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_exp_outerlayers         ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_sta_lostHits                ( new vector<int>           );
    unique_ptr<vector<float> >         vector_mus_sta_ptErr                   ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_etaErr                  ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_sta_phiErr                  ( new vector<float>         );
    unique_ptr<vector<int> >           vector_mus_sta_trk_charge              ( new vector<int>           );

    ////////////  
    // Global //
    ////////////
    unique_ptr<vector<float> >         vector_mus_gfit_d0                     ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_gfit_z0                     ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_gfit_d0Err                  ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_gfit_z0Err                  ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_gfit_d0corr                 ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_gfit_z0corr                 ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_gfit_qoverp                 ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_gfit_qoverpError            ( new vector<float>         );
    unique_ptr<vector<int> >           vector_mus_gfit_validHits              ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_gfit_validSiHits            ( new vector<int>           );
    unique_ptr<vector<LorentzVector> > vector_mus_gfit_vertex_p4              ( new vector<LorentzVector> );
    unique_ptr<vector<float> >         vector_mus_gfit_d0corrPhi              ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_gfit_d0phiCov               ( new vector<float>         );
    unique_ptr<vector<int> >           vector_mus_gfit_qualityMask            ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_gfit_nlayers                ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_gfit_nlayers3D              ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_gfit_nlayersLost            ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_gfit_validPixelHits         ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_gfit_lostPixelHits          ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_gfit_exp_innerlayers        ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_gfit_exp_outerlayers        ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_gfit_lostHits               ( new vector<int>           );
    unique_ptr<vector<float> >         vector_mus_gfit_etaErr                 ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_gfit_phiErr                 ( new vector<float>         );
    unique_ptr<vector<int> >           vector_mus_gfit_trk_charge             ( new vector<int>           );   

    ////////////  
    // Best   //
    ////////////

    unique_ptr<vector<float> >         vector_mus_bfit_d0                     ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_z0                     ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_d0Err                  ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_z0Err                  ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_d0corr                 ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_z0corr                 ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_qoverp                 ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_qoverpError            ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_chi2                   ( new vector<float>         );
    unique_ptr<vector<int> >           vector_mus_bfit_ndof                   ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_validHits              ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_validSTAHits           ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_validSiHits            ( new vector<int>           );
    unique_ptr<vector<LorentzVector> > vector_mus_bfit_vertex_p4              ( new vector<LorentzVector> );
    unique_ptr<vector<float> >         vector_mus_bfit_d0corrPhi              ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_d0phiCov               ( new vector<float>         );
    unique_ptr<vector<int> >           vector_mus_bfit_qualityMask            ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_nlayers                ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_nlayers3D              ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_nlayersLost            ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_validPixelHits         ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_lostPixelHits          ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_exp_innerlayers        ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_exp_outerlayers        ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_bfit_lostHits               ( new vector<int>           );
    unique_ptr<vector<float> >         vector_mus_bfit_etaErr                 ( new vector<float>         );
    unique_ptr<vector<float> >         vector_mus_bfit_phiErr                 ( new vector<float>         );
    unique_ptr<vector<int> >           vector_mus_bfit_trk_charge             ( new vector<int>           );   

    /////////////
    // Quality //
    /////////////

    unique_ptr<vector<bool> >   vector_mus_updatedSta          ( new vector<bool>   );
    unique_ptr<vector<bool> >   vector_mus_tightMatch          ( new vector<bool>   );
    unique_ptr<vector<float> >  vector_mus_glbKink             ( new vector<float>  );
    unique_ptr<vector<float> >  vector_mus_trkRelChi2          ( new vector<float>  );
    unique_ptr<vector<float> >  vector_mus_staRelChi2          ( new vector<float>  );
    unique_ptr<vector<float> >  vector_mus_localDistance       ( new vector<float>  );
    unique_ptr<vector<float> >  vector_mus_globalDeltaEtaPhi   ( new vector<float>  );
    unique_ptr<vector<float> >  vector_mus_glbTrackProbability ( new vector<float>  );


    ////////////
    // Timing //
    ////////////
  
    unique_ptr<vector<int> >           vector_mus_timeNumStationsUsed         ( new vector<int>     );
    unique_ptr<vector<float> >         vector_mus_timeAtIpInOut               ( new vector<float>   );
    unique_ptr<vector<float> >         vector_mus_timeAtIpInOutErr            ( new vector<float>   );
    unique_ptr<vector<float> >         vector_mus_timeAtIpOutIn               ( new vector<float>   );
    unique_ptr<vector<float> >         vector_mus_timeAtIpOutInErr            ( new vector<float>   );
    unique_ptr<vector<int> >           vector_mus_timeDirection               ( new vector<int>     );


    ///////////
    // Muons //
    ///////////

    unique_ptr<vector<int> >           vector_mus_goodmask                ( new vector<int>           );        
    unique_ptr<vector<int> >           vector_mus_nmatches                ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_charge                  ( new vector<int>           );        
    unique_ptr<vector<int> >           vector_mus_nOverlaps               ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_overlap0                ( new vector<int>           );
    unique_ptr<vector<int> >           vector_mus_overlap1                ( new vector<int>           );
    unique_ptr<vector<float> >         vector_mus_mass                    ( new vector<float>         );
    unique_ptr<vector<bool> >          vector_mus_isRPCMuon               ( new vector<bool>          );


    ////////////
    // Energy //
    ////////////

    unique_ptr<vector<float> >         vector_mus_e_em                ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_e_had               ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_e_ho                ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_e_emS9              ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_e_hadS9             ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_e_hoS9              ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_emS25               ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_emMax               ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_hadMax              ( new vector<float>          );
    unique_ptr<vector<int> >           vector_mus_ecal_rawId          ( new vector<int>            );
    unique_ptr<vector<int> >           vector_mus_hcal_rawId          ( new vector<int>            );


    ///////////////
    // Isolation //
    ///////////////

    unique_ptr<vector<float> >         vector_mus_iso03_hoEt          ( new vector<float>          );
    unique_ptr<vector<int> >           vector_mus_iso03_njets         ( new vector<int>            );
    unique_ptr<vector<float> >         vector_mus_iso05_sumPt         ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_iso05_emEt          ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_iso05_hadEt         ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_iso05_hoEt          ( new vector<float>          );
    unique_ptr<vector<int> >           vector_mus_iso05_ntrk          ( new vector<int>            );


    ////////////
    // Tracks //
    ////////////

    unique_ptr<vector<int> >           vector_mus_muonBestTrackType   ( new vector<int>            );        
    unique_ptr<vector<LorentzVector> > vector_mus_vertex_p4           ( new vector<LorentzVector>  );
    unique_ptr<vector<float> >         vector_mus_d0                  ( new vector<float>          );      
    unique_ptr<vector<float> >         vector_mus_z0                  ( new vector<float>          );      
    unique_ptr<vector<float> >         vector_mus_d0corr              ( new vector<float>          );      
    unique_ptr<vector<float> >         vector_mus_z0corr              ( new vector<float>          );      
    unique_ptr<vector<float> >         vector_mus_vertexphi           ( new vector<float>          );      
    unique_ptr<vector<float> >         vector_mus_chi2                ( new vector<float>          );      
    unique_ptr<vector<float> >         vector_mus_ndof                ( new vector<float>          );      
    unique_ptr<vector<float> >         vector_mus_etaErr              ( new vector<float>          );      
    unique_ptr<vector<float> >         vector_mus_phiErr              ( new vector<float>          );      
    unique_ptr<vector<int> >           vector_mus_trk_charge          ( new vector<int>            );   
    unique_ptr<vector<float> >         vector_mus_qoverp              ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_qoverpError         ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_d0corrPhi           ( new vector<float>          );
    unique_ptr<vector<float> >         vector_mus_d0phiCov            ( new vector<float>          );
    unique_ptr<vector<int> >           vector_mus_qualityMask         ( new vector<int>            );
    unique_ptr<vector<int> >           vector_mus_nlayers3D           ( new vector<int>            );
    unique_ptr<vector<int> >           vector_mus_nlayersLost         ( new vector<int>            );
    unique_ptr<vector<int> >           vector_mus_lostPixelHits       ( new vector<int>            );

  
    ////////
    // PF //
    ////////

    unique_ptr< vector<float> >         vector_mus_isoMeanDRR03_pf_ChargedHadronPt                 ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR03_pf_ChargedParticlePt               ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR03_pf_NeutralHadronEt                 ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR03_pf_PhotonEt                        ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR03_pf_sumNeutralHadronEtHighThreshold ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR03_pf_sumPhotonEtHighThreshold        ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR03_pf_PUPt                            ( new vector<float>   );

    unique_ptr< vector<float> >         vector_mus_isoMeanDRR04_pf_ChargedHadronPt                 ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR04_pf_ChargedParticlePt               ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR04_pf_NeutralHadronEt                 ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR04_pf_PhotonEt                        ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR04_pf_sumNeutralHadronEtHighThreshold ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR04_pf_sumPhotonEtHighThreshold        ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoMeanDRR04_pf_PUPt                            ( new vector<float>   );

    unique_ptr< vector<float> >         vector_mus_isoSumDRR03_pf_ChargedHadronPt                 ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR03_pf_ChargedParticlePt               ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR03_pf_NeutralHadronEt                 ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR03_pf_PhotonEt                        ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR03_pf_sumNeutralHadronEtHighThreshold ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR03_pf_sumPhotonEtHighThreshold        ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR03_pf_PUPt                            ( new vector<float>   );

    unique_ptr< vector<float> >         vector_mus_isoSumDRR04_pf_ChargedHadronPt                 ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR04_pf_ChargedParticlePt               ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR04_pf_NeutralHadronEt                 ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR04_pf_PhotonEt                        ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR04_pf_sumNeutralHadronEtHighThreshold ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR04_pf_sumPhotonEtHighThreshold        ( new vector<float>   );
    unique_ptr< vector<float> >         vector_mus_isoSumDRR04_pf_PUPt                            ( new vector<float>   );

    ///////////
    // IP 3D //
    ///////////

    unique_ptr<vector<float> >         vector_mus_bs3d                        ( new vector<float>   );
    unique_ptr<vector<float> >         vector_mus_bs3derr                     ( new vector<float>   );
    unique_ptr<vector<float> >         vector_mus_bs2d                        ( new vector<float>   );
    unique_ptr<vector<float> >         vector_mus_bs2derr                     ( new vector<float>   );


    ////////////////////////////
    // --- Fill Muon Data --- //
    ////////////////////////////

    ///////////////
    // Get Muons //
    ///////////////

    Handle<View<pat::Muon> > muon_h;
    iEvent.getByToken( muonsToken , muon_h );

    //////////////////
    // Get Vertices //
    //////////////////

    iEvent.getByToken( vtxToken , vertexHandle );

    ///////////////////////
    // Get PF Candidates //
    ///////////////////////

    iEvent.getByToken(pfCandsToken, packPfCand_h);
    pfCandidates  = packPfCand_h.product();


    /////////////////////////////////////
    // Get BeamSpot from BeamSpotMaker //
    /////////////////////////////////////

    Handle<LorentzVector> beamSpotH;
    iEvent.getByToken(beamSpotToken, beamSpotH);
    const Point beamSpot = beamSpotH.isValid() ? Point(beamSpotH->x(), beamSpotH->y(), beamSpotH->z()) : Point(0,0,0);
  
    ///////////
    // Muons // 
    ///////////
  
    unsigned int muonIndex = 0;
    View<pat::Muon>::const_iterator muons_end = muon_h->end();  // Iterator
    for ( View<pat::Muon>::const_iterator muon = muon_h->begin(); muon != muons_end; ++muon ) {


        // References
        const RefToBase<pat::Muon>    muonRef                 = muon_h->refAt(muonIndex); 
        const TrackRef                globalTrack             = muon->globalTrack();
        const TrackRef                siTrack                 = muon->innerTrack();
        const TrackRef                staTrack                = muon->outerTrack();
        const TrackRef                bestTrack               = muon->muonBestTrack();
        const MuonQuality             quality                 = muon->combinedQuality();
        const VertexCollection*       vertexCollection        = vertexHandle.product();

        // Iterators
        VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
        int firstGoodVertexIdx = 0;
        for (VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++firstGoodVertexIdx) {
            if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0 ) {
                firstGoodVertex = vtx;
                break;
            }
        }

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
        vector_mus_sta_p4            -> push_back( staTrack.isNonnull()  ? LorentzVector( staTrack->px(), staTrack->py(), staTrack->pz(), staTrack->p() ) : LorentzVector(0, 0, 0, 0) );
        vector_mus_sta_vertex_p4     -> push_back( staTrack.isNonnull()  ? LorentzVector( staTrack->vx(), staTrack->vy(), staTrack->vz(), 0.0           ) : LorentzVector(0, 0, 0, 0) );
        // Embedding all trackMaker details
        vector_mus_sta_d0corrPhi       -> push_back( staTrack.isNonnull() ? atan2( (    staTrack->dxy(beamSpot) * sin( staTrack->phi() ) ),
                                                                                   -1 * staTrack->dxy(beamSpot) * cos( staTrack->phi() ) )                       : -9999. );
        vector_mus_sta_d0phiCov        -> push_back( staTrack.isNonnull() ? -1.* staTrack->covariance(TrackBase::i_phi, TrackBase::i_dxy)                        : -9999. );
        vector_mus_sta_qualityMask     -> push_back( staTrack.isNonnull() ? staTrack->qualityMask()                                                              : -9999  );
        vector_mus_sta_algo            -> push_back( staTrack.isNonnull() ? staTrack->algo()                                                                     : -9999  );
        vector_mus_sta_nlayers         -> push_back( staTrack.isNonnull() ? staTrack->hitPattern().trackerLayersWithMeasurement()                                : -9999  );
        vector_mus_sta_nlayers3D       -> push_back( staTrack.isNonnull() ? (staTrack->hitPattern().pixelLayersWithMeasurement() +
                                                                             staTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo())                 : -9999  );
        vector_mus_sta_nlayersLost     -> push_back( staTrack.isNonnull() ? staTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) : -9999  );
        vector_mus_sta_validPixelHits  -> push_back( staTrack.isNonnull() ? staTrack->hitPattern().numberOfValidPixelHits()                                      : -9999  );
        vector_mus_sta_lostPixelHits   -> push_back( staTrack.isNonnull() ? staTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS)           : -9999  );
        vector_mus_sta_exp_innerlayers -> push_back( staTrack.isNonnull() ? staTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS)            : -9999  );
        vector_mus_sta_exp_outerlayers -> push_back( staTrack.isNonnull() ? staTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS)            : -9999  );
        vector_mus_sta_lostHits        -> push_back( staTrack.isNonnull() ? staTrack->numberOfLostHits()                                                         : -9999  );
        vector_mus_sta_ptErr           -> push_back( staTrack.isNonnull() ? staTrack->ptError()                                                                  : -9999. );
        vector_mus_sta_etaErr          -> push_back( staTrack.isNonnull() ? staTrack->etaError()                                                                 : -9999. );
        vector_mus_sta_phiErr          -> push_back( staTrack.isNonnull() ? staTrack->phiError()                                                                 : -9999. );
        vector_mus_sta_trk_charge      -> push_back( staTrack.isNonnull() ? staTrack->charge()                                                                   : -9999  );

        ////////////
        // Global //
        ////////////
        vector_mus_gfit_d0           -> push_back( globalTrack.isNonnull() ? globalTrack->d0()                                    : -9999. );
        vector_mus_gfit_z0           -> push_back( globalTrack.isNonnull() ? globalTrack->dz()                                    : -9999. );
        vector_mus_gfit_d0Err        -> push_back( globalTrack.isNonnull() ? globalTrack->d0Error()                               : -9999. );
        vector_mus_gfit_z0Err        -> push_back( globalTrack.isNonnull() ? globalTrack->dzError()                               : -9999. );
        vector_mus_gfit_d0corr       -> push_back( globalTrack.isNonnull() ? -1*(globalTrack->dxy(beamSpot))                      : -9999. );
        vector_mus_gfit_z0corr       -> push_back( globalTrack.isNonnull() ? globalTrack->dz(beamSpot)                            : -9999. );
        vector_mus_gfit_qoverp       -> push_back( globalTrack.isNonnull() ? globalTrack->qoverp()                                : -9999. );
        vector_mus_gfit_qoverpError  -> push_back( globalTrack.isNonnull() ? globalTrack->qoverpError()                           : -9999. );
        vector_mus_gfit_validHits    -> push_back( globalTrack.isNonnull() ? globalTrack->numberOfValidHits()                     : -9999  );
        vector_mus_gfit_validSiHits  -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfValidTrackerHits() : -9999  );
        vector_mus_gfit_vertex_p4    -> push_back( globalTrack.isNonnull() ? LorentzVector( globalTrack->vx(), globalTrack->vy(), globalTrack->vz(), 0.0              ) : LorentzVector(0.0,0.0,0.0,0.0) );
	// Embedding all trackMaker details
        vector_mus_gfit_d0corrPhi       -> push_back( globalTrack.isNonnull() ? atan2( (    globalTrack->dxy(beamSpot) * sin( globalTrack->phi() ) ),
                                                                                       -1 * globalTrack->dxy(beamSpot) * cos( globalTrack->phi() ) )                    : -9999. );
        vector_mus_gfit_d0phiCov        -> push_back( globalTrack.isNonnull() ? -1.* globalTrack->covariance(TrackBase::i_phi, TrackBase::i_dxy)                        : -9999. );
        vector_mus_gfit_qualityMask     -> push_back( globalTrack.isNonnull() ? globalTrack->qualityMask()                                                              : -9999  );
        vector_mus_gfit_nlayers         -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().trackerLayersWithMeasurement()                                : -9999  );
        vector_mus_gfit_nlayers3D       -> push_back( globalTrack.isNonnull() ? (globalTrack->hitPattern().pixelLayersWithMeasurement() +
                                                                                 globalTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo())                 : -9999  );
        vector_mus_gfit_nlayersLost     -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) : -9999  );
        vector_mus_gfit_validPixelHits  -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfValidPixelHits()                                      : -9999  );
        vector_mus_gfit_lostPixelHits   -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS)           : -9999  );
        vector_mus_gfit_exp_innerlayers -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS)            : -9999  );
        vector_mus_gfit_exp_outerlayers -> push_back( globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS)            : -9999  );
        vector_mus_gfit_lostHits        -> push_back( globalTrack.isNonnull() ? globalTrack->numberOfLostHits()                                                         : -9999  );
        vector_mus_gfit_etaErr          -> push_back( globalTrack.isNonnull() ? globalTrack->etaError()                                                                 : -9999. );
        vector_mus_gfit_phiErr          -> push_back( globalTrack.isNonnull() ? globalTrack->phiError()                                                                 : -9999. );
        vector_mus_gfit_trk_charge      -> push_back( globalTrack.isNonnull() ? globalTrack->charge()                                                                   : -9999  );


	//////////
        // Best //
        //////////
        vector_mus_bfit_d0           -> push_back( bestTrack.isNonnull() ? bestTrack->d0()                                    : -9999. );
        vector_mus_bfit_z0           -> push_back( bestTrack.isNonnull() ? bestTrack->dz()                                    : -9999. );
        vector_mus_bfit_d0Err        -> push_back( bestTrack.isNonnull() ? bestTrack->d0Error()                               : -9999. );
        vector_mus_bfit_z0Err        -> push_back( bestTrack.isNonnull() ? bestTrack->dzError()                               : -9999. );
        vector_mus_bfit_d0corr       -> push_back( bestTrack.isNonnull() ? -1*(bestTrack->dxy(beamSpot))                      : -9999. );
        vector_mus_bfit_z0corr       -> push_back( bestTrack.isNonnull() ? bestTrack->dz(beamSpot)                            : -9999. );
        vector_mus_bfit_qoverp       -> push_back( bestTrack.isNonnull() ? bestTrack->qoverp()                                : -9999. );
        vector_mus_bfit_qoverpError  -> push_back( bestTrack.isNonnull() ? bestTrack->qoverpError()                           : -9999. );
        vector_mus_bfit_chi2         -> push_back( bestTrack.isNonnull() ? bestTrack->chi2()                                  : -9999. );
        vector_mus_bfit_ndof         -> push_back( bestTrack.isNonnull() ? bestTrack->ndof()                                  : -9999  );
        vector_mus_bfit_validHits    -> push_back( bestTrack.isNonnull() ? bestTrack->numberOfValidHits()                     : -9999  );
        vector_mus_bfit_validSTAHits -> push_back( bestTrack.isNonnull() ? bestTrack->hitPattern().numberOfValidMuonHits()    : -9999  );
        vector_mus_bfit_validSiHits  -> push_back( bestTrack.isNonnull() ? bestTrack->hitPattern().numberOfValidTrackerHits() : -9999  );
        vector_mus_bfit_vertex_p4    -> push_back( bestTrack.isNonnull() ? LorentzVector( bestTrack->vx(), bestTrack->vy(), bestTrack->vz(),  0.0           ) : LorentzVector(0.0,0.0,0.0,0.0) );
	// Embedding all trackMaker details
        vector_mus_bfit_d0corrPhi       -> push_back( bestTrack.isNonnull() ? atan2( (    bestTrack->dxy(beamSpot) * sin( bestTrack->phi() ) ),
                                                                                     -1 * bestTrack->dxy(beamSpot) * cos( bestTrack->phi() ) )                      : -9999. );
        vector_mus_bfit_d0phiCov        -> push_back( bestTrack.isNonnull() ? -1.* bestTrack->covariance(TrackBase::i_phi, TrackBase::i_dxy)                        : -9999. );
        vector_mus_bfit_qualityMask     -> push_back( bestTrack.isNonnull() ? bestTrack->qualityMask()                                                              : -9999. );
        vector_mus_bfit_nlayers         -> push_back( bestTrack.isNonnull() ? bestTrack->hitPattern().trackerLayersWithMeasurement()                                : -9999  );
        vector_mus_bfit_nlayers3D       -> push_back( bestTrack.isNonnull() ? (bestTrack->hitPattern().pixelLayersWithMeasurement() +
                                                                               bestTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo())                 : -9999  );
        vector_mus_bfit_nlayersLost     -> push_back( bestTrack.isNonnull() ? bestTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) : -9999  );
        vector_mus_bfit_validPixelHits  -> push_back( bestTrack.isNonnull() ? bestTrack->hitPattern().numberOfValidPixelHits()                                      : -9999  );
        vector_mus_bfit_lostPixelHits   -> push_back( bestTrack.isNonnull() ? bestTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS)           : -9999  );
        vector_mus_bfit_exp_innerlayers -> push_back( bestTrack.isNonnull() ? bestTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS)            : -9999  );
        vector_mus_bfit_exp_outerlayers -> push_back( bestTrack.isNonnull() ? bestTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS)            : -9999  );
        vector_mus_bfit_lostHits        -> push_back( bestTrack.isNonnull() ? bestTrack->numberOfLostHits()                                                         : -9999  );
        vector_mus_bfit_etaErr          -> push_back( bestTrack.isNonnull() ? bestTrack->etaError()                                                                 : -9999. );
        vector_mus_bfit_phiErr          -> push_back( bestTrack.isNonnull() ? bestTrack->phiError()                                                                 : -9999. );
        vector_mus_bfit_trk_charge      -> push_back( bestTrack.isNonnull() ? bestTrack->charge()                                                                   : -9999  );

        //////////////////
        // Muon Quality //
        //////////////////
        vector_mus_updatedSta          -> push_back( quality.updatedSta          );
        vector_mus_tightMatch          -> push_back( quality.tightMatch          );
        vector_mus_glbKink             -> push_back( quality.glbKink             );
        vector_mus_trkRelChi2          -> push_back( quality.trkRelChi2          );
        vector_mus_staRelChi2          -> push_back( quality.staRelChi2          );
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

        //////////
        // Muon //
        //////////

        // Calculate Overlaps
        int mus_overlap0 = -1, mus_overlap1 = -1, muInd = -1, mus_nOverlaps = 0;
        for ( View<pat::Muon>::const_iterator muonJ = muon_h->begin(); muonJ != muons_end; ++muonJ ) {
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

        vector_mus_goodmask                -> push_back( goodMask                                                  );
        vector_mus_nmatches                -> push_back( muon->isMatchesValid() ? muon->numberOfMatches() :  -9999 );
        vector_mus_nOverlaps               -> push_back( mus_nOverlaps                                             );
        vector_mus_overlap0                -> push_back( mus_overlap0                                              );
        vector_mus_overlap1                -> push_back( mus_overlap1                                              );
        vector_mus_mass                    -> push_back( muon->mass()                                              );
        vector_mus_isRPCMuon               ->push_back( muon->isRPCMuon()                                          );


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

        //calEnergy is an object of type reco::MuonEnergy, learn about it at
        //cmslxr.fnal.gov/lxr/source/DataFormats/MuonReco/interface/MuonEnergy.h
        vector_mus_emS25           -> push_back( energyIsValid ? muon->calEnergy().emS25                         : -9999. );
        vector_mus_emMax           -> push_back( energyIsValid ? muon->calEnergy().emMax                         : -9999. );
        vector_mus_hadMax          -> push_back( energyIsValid ? muon->calEnergy().hadMax                        : -9999. );
        vector_mus_ecal_rawId      -> push_back( energyIsValid ? muon->calEnergy().ecal_id.rawId()               : -9999. );
        vector_mus_hcal_rawId      -> push_back( energyIsValid ? muon->calEnergy().hcal_id.rawId()               : -9999. );
    

        ///////////////
        // Isolation //
        ///////////////

        vector_mus_iso03_hoEt         -> push_back( muon->isIsolationValid() ? muon->isolationR03().hoEt           : -9999.        );
        vector_mus_iso03_njets        -> push_back( muon->isIsolationValid() ? muon->isolationR03().nJets          : -9999         );
        vector_mus_iso05_sumPt        -> push_back( muon->isIsolationValid() ? muon->isolationR05().sumPt          : -9999.        );
        vector_mus_iso05_emEt         -> push_back( muon->isIsolationValid() ? muon->isolationR05().emEt           : -9999.        );
        vector_mus_iso05_hadEt        -> push_back( muon->isIsolationValid() ? muon->isolationR05().hadEt          : -9999.        );
        vector_mus_iso05_hoEt         -> push_back( muon->isIsolationValid() ? muon->isolationR05().hoEt           : -9999.        );
        vector_mus_iso05_ntrk         -> push_back( muon->isIsolationValid() ? muon->isolationR05().nTracks        : -9999         );


        ////////////
        // Tracks //
        ////////////

        vector_mus_muonBestTrackType  -> push_back( muon->muonBestTrackType() );
        vector_mus_vertex_p4          -> push_back( siTrack.isNonnull()     ? LorentzVector( siTrack->vx()       , siTrack->vy()       , siTrack->vz()       , 0.0                ) : LorentzVector( -9999.0, -9999.0, -9999.0, -9999.0) );
        vector_mus_d0                 -> push_back( siTrack.isNonnull()     ? siTrack->d0()                                        : -9999.        );
        vector_mus_z0                 -> push_back( siTrack.isNonnull()     ? siTrack->dz()                                        : -9999.        );
        vector_mus_d0corr             -> push_back( siTrack.isNonnull()     ? -1*(siTrack->dxy(beamSpot))                          : -9999.        );
        vector_mus_z0corr             -> push_back( siTrack.isNonnull()     ? siTrack->dz(beamSpot)                                : -9999.        );
        vector_mus_vertexphi          -> push_back( siTrack.isNonnull()     ? atan2( siTrack->vy(), siTrack->vx() )                : -9999.        );
        vector_mus_chi2               -> push_back( siTrack.isNonnull()     ? siTrack->chi2()                                      : -9999.        );
        vector_mus_ndof               -> push_back( siTrack.isNonnull()     ? siTrack->ndof()                                      : -9999.        );
        vector_mus_etaErr             -> push_back( siTrack.isNonnull()     ? siTrack->etaError()                                  :  -9999.       );
        vector_mus_phiErr             -> push_back( siTrack.isNonnull()     ? siTrack->phiError()                                  :  -9999.       );
        vector_mus_trk_charge         -> push_back( siTrack.isNonnull()     ? siTrack->charge()                                    :  -9999        );
        vector_mus_qoverp             -> push_back( siTrack.isNonnull()     ? siTrack->qoverp()                                    :  -9999.       );
        vector_mus_qoverpError        -> push_back( siTrack.isNonnull()     ? siTrack->qoverpError()                               :  -9999.       );
        // Embedding all trackMaker details
        vector_mus_d0corrPhi          -> push_back( siTrack.isNonnull()     ?  atan2( (siTrack->dxy(beamSpot) * sin( siTrack->phi() )), -1 * siTrack->dxy(beamSpot) * cos( siTrack->phi() ) ) : -9999.        );
        vector_mus_d0phiCov           -> push_back( siTrack.isNonnull()     ?  -1.* siTrack->covariance(TrackBase::i_phi, TrackBase::i_dxy) : -9999. );
        vector_mus_qualityMask        -> push_back( siTrack.isNonnull()     ? siTrack->qualityMask()                               : -9999.        );
        vector_mus_nlayers3D          -> push_back( siTrack.isNonnull()     ? siTrack->hitPattern().pixelLayersWithMeasurement()  + siTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo():  -9999        );
        vector_mus_nlayersLost        -> push_back( siTrack.isNonnull()     ? siTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) :  -9999     );
        vector_mus_lostPixelHits      -> push_back( siTrack.isNonnull()     ? siTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS)        :  -9999        );

        ////////
        // PF //
        ////////

        // PF Isolation
        MuonPFIsolation pfStructMeanDRR03 = muon->pfMeanDRIsoProfileR03();
        MuonPFIsolation pfStructMeanDRR04 = muon->pfSumDRIsoProfileR04();
        MuonPFIsolation pfStructSumDRR03  = muon->pfMeanDRIsoProfileR03();
        MuonPFIsolation pfStructSumDRR04  = muon->pfSumDRIsoProfileR04();

        vector_mus_isoMeanDRR03_pf_ChargedHadronPt                 -> push_back( pfStructMeanDRR03.sumChargedHadronPt              );
        vector_mus_isoMeanDRR03_pf_ChargedParticlePt               -> push_back( pfStructMeanDRR03.sumChargedParticlePt            );
        vector_mus_isoMeanDRR03_pf_NeutralHadronEt                 -> push_back( pfStructMeanDRR03.sumNeutralHadronEt              );
        vector_mus_isoMeanDRR03_pf_PhotonEt                        -> push_back( pfStructMeanDRR03.sumPhotonEt                     );
        vector_mus_isoMeanDRR03_pf_sumNeutralHadronEtHighThreshold -> push_back( pfStructMeanDRR03.sumNeutralHadronEtHighThreshold );
        vector_mus_isoMeanDRR03_pf_sumPhotonEtHighThreshold        -> push_back( pfStructMeanDRR03.sumPhotonEtHighThreshold        );
        vector_mus_isoMeanDRR03_pf_PUPt                            -> push_back( pfStructMeanDRR03.sumPUPt                         );

        vector_mus_isoMeanDRR04_pf_ChargedHadronPt                 -> push_back( pfStructMeanDRR04.sumChargedHadronPt              );
        vector_mus_isoMeanDRR04_pf_ChargedParticlePt               -> push_back( pfStructMeanDRR04.sumChargedParticlePt            );
        vector_mus_isoMeanDRR04_pf_NeutralHadronEt                 -> push_back( pfStructMeanDRR04.sumNeutralHadronEt              );
        vector_mus_isoMeanDRR04_pf_PhotonEt                        -> push_back( pfStructMeanDRR04.sumPhotonEt                     );
        vector_mus_isoMeanDRR04_pf_sumNeutralHadronEtHighThreshold -> push_back( pfStructMeanDRR04.sumNeutralHadronEtHighThreshold );
        vector_mus_isoMeanDRR04_pf_sumPhotonEtHighThreshold        -> push_back( pfStructMeanDRR04.sumPhotonEtHighThreshold        );
        vector_mus_isoMeanDRR04_pf_PUPt                            -> push_back( pfStructMeanDRR04.sumPUPt                         );

        vector_mus_isoSumDRR03_pf_ChargedHadronPt                 -> push_back( pfStructSumDRR03.sumChargedHadronPt              );
        vector_mus_isoSumDRR03_pf_ChargedParticlePt               -> push_back( pfStructSumDRR03.sumChargedParticlePt            );
        vector_mus_isoSumDRR03_pf_NeutralHadronEt                 -> push_back( pfStructSumDRR03.sumNeutralHadronEt              );
        vector_mus_isoSumDRR03_pf_PhotonEt                        -> push_back( pfStructSumDRR03.sumPhotonEt                     );
        vector_mus_isoSumDRR03_pf_sumNeutralHadronEtHighThreshold -> push_back( pfStructSumDRR03.sumNeutralHadronEtHighThreshold );
        vector_mus_isoSumDRR03_pf_sumPhotonEtHighThreshold        -> push_back( pfStructSumDRR03.sumPhotonEtHighThreshold        );
        vector_mus_isoSumDRR03_pf_PUPt                            -> push_back( pfStructSumDRR03.sumPUPt                         );

        vector_mus_isoSumDRR04_pf_ChargedHadronPt                 -> push_back( pfStructSumDRR04.sumChargedHadronPt              );
        vector_mus_isoSumDRR04_pf_ChargedParticlePt               -> push_back( pfStructSumDRR04.sumChargedParticlePt            );
        vector_mus_isoSumDRR04_pf_NeutralHadronEt                 -> push_back( pfStructSumDRR04.sumNeutralHadronEt              );
        vector_mus_isoSumDRR04_pf_PhotonEt                        -> push_back( pfStructSumDRR04.sumPhotonEt                     );
        vector_mus_isoSumDRR04_pf_sumNeutralHadronEtHighThreshold -> push_back( pfStructSumDRR04.sumNeutralHadronEtHighThreshold );
        vector_mus_isoSumDRR04_pf_sumPhotonEtHighThreshold        -> push_back( pfStructSumDRR04.sumPhotonEtHighThreshold        );
        vector_mus_isoSumDRR04_pf_PUPt                            -> push_back( pfStructSumDRR04.sumPUPt                         );

        ///////////
        // IP 3D //
        ///////////

        vector_mus_bs3d         -> push_back( muon->dB(pat::Muon::BS3D) ); 
        vector_mus_bs3derr      -> push_back( muon->edB(pat::Muon::BS3D) );
        vector_mus_bs2d         -> push_back( muon->dB(pat::Muon::BS2D) ); 
        vector_mus_bs2derr      -> push_back( muon->edB(pat::Muon::BS2D) );
    
        muonIndex++;

    } // end loop on muons

    /////////
    // STA //
    /////////

    iEvent.put(std::move( vector_mus_sta_d0                       ), branchprefix_ + "stad0"              );
    iEvent.put(std::move( vector_mus_sta_z0                       ), branchprefix_ + "staz0"              );
    iEvent.put(std::move( vector_mus_sta_d0Err                    ), branchprefix_ + "stad0Err"           );
    iEvent.put(std::move( vector_mus_sta_z0Err                    ), branchprefix_ + "staz0Err"           );
    iEvent.put(std::move( vector_mus_sta_d0corr                   ), branchprefix_ + "stad0corr"          );
    iEvent.put(std::move( vector_mus_sta_z0corr                   ), branchprefix_ + "staz0corr"          );
    iEvent.put(std::move( vector_mus_sta_qoverp                   ), branchprefix_ + "staqoverp"          );
    iEvent.put(std::move( vector_mus_sta_qoverpError              ), branchprefix_ + "staqoverpError"     );
    iEvent.put(std::move( vector_mus_sta_chi2                     ), branchprefix_ + "stachi2"            );
    iEvent.put(std::move( vector_mus_sta_ndof                     ), branchprefix_ + "standof"            );
    iEvent.put(std::move( vector_mus_sta_validHits                ), branchprefix_ + "stavalidHits"       );
    iEvent.put(std::move( vector_mus_sta_p4                       ), branchprefix_ + "stap4"              );
    iEvent.put(std::move( vector_mus_sta_vertex_p4                ), branchprefix_ + "stavertexp4"        );
    iEvent.put(std::move( vector_mus_sta_d0corrPhi                ), branchprefix_ + "stad0corrPhi"       );
    iEvent.put(std::move( vector_mus_sta_d0phiCov        	      ), branchprefix_ + "stad0phiCov"        );
    iEvent.put(std::move( vector_mus_sta_qualityMask     	      ), branchprefix_ + "staqualityMask"     );
    iEvent.put(std::move( vector_mus_sta_algo              	      ), branchprefix_ + "staalgo"            );
    iEvent.put(std::move( vector_mus_sta_nlayers         	      ), branchprefix_ + "stanlayers"         );
    iEvent.put(std::move( vector_mus_sta_nlayers3D       	      ), branchprefix_ + "stanlayers3D"       );
    iEvent.put(std::move( vector_mus_sta_nlayersLost     	      ), branchprefix_ + "stanlayersLost"     );
    iEvent.put(std::move( vector_mus_sta_validPixelHits  	      ), branchprefix_ + "stavalidPixelHits"  );
    iEvent.put(std::move( vector_mus_sta_lostPixelHits   	      ), branchprefix_ + "stalostPixelHits"   );
    iEvent.put(std::move( vector_mus_sta_exp_innerlayers 	      ), branchprefix_ + "staexpinnerlayer"  );
    iEvent.put(std::move( vector_mus_sta_exp_outerlayers 	      ), branchprefix_ + "staexpouterlayer"  );
    iEvent.put(std::move( vector_mus_sta_lostHits        	      ), branchprefix_ + "stalostHits"        );
    iEvent.put(std::move( vector_mus_sta_ptErr           	      ), branchprefix_ + "staptErr"           );
    iEvent.put(std::move( vector_mus_sta_etaErr          	      ), branchprefix_ + "staetaErr"          );
    iEvent.put(std::move( vector_mus_sta_phiErr          	      ), branchprefix_ + "staphiErr"          );
    iEvent.put(std::move( vector_mus_sta_trk_charge               ), branchprefix_ + "statrkcharge"       );

    ////////////                                                                       
    // Global //
    /////////////

    iEvent.put(std::move( vector_mus_gfit_d0                      ), branchprefix_ + "gfitd0"             );
    iEvent.put(std::move( vector_mus_gfit_z0                      ), branchprefix_ + "gfitz0"             );
    iEvent.put(std::move( vector_mus_gfit_d0Err                   ), branchprefix_ + "gfitd0Err"          );
    iEvent.put(std::move( vector_mus_gfit_z0Err                   ), branchprefix_ + "gfitz0Err"          );
    iEvent.put(std::move( vector_mus_gfit_d0corr                  ), branchprefix_ + "gfitd0corr"         );
    iEvent.put(std::move( vector_mus_gfit_z0corr                  ), branchprefix_ + "gfitz0corr"         );
    iEvent.put(std::move( vector_mus_gfit_qoverp                  ), branchprefix_ + "gfitqoverp"         );
    iEvent.put(std::move( vector_mus_gfit_qoverpError             ), branchprefix_ + "gfitqoverpError"    );
    iEvent.put(std::move( vector_mus_gfit_validHits               ), branchprefix_ + "gfitvalidHits"      );
    iEvent.put(std::move( vector_mus_gfit_validSiHits             ), branchprefix_ + "gfitvalidSiHits"    );
    iEvent.put(std::move( vector_mus_gfit_vertex_p4               ), branchprefix_ + "gfitvertexp4"             );
    iEvent.put(std::move( vector_mus_gfit_d0corrPhi               ), branchprefix_ + "gfitd0corrPhi"       );
    iEvent.put(std::move( vector_mus_gfit_d0phiCov        	      ), branchprefix_ + "gfitd0phiCov"        );
    iEvent.put(std::move( vector_mus_gfit_qualityMask     	      ), branchprefix_ + "gfitqualityMask"     );
    iEvent.put(std::move( vector_mus_gfit_nlayers         	      ), branchprefix_ + "gfitnlayers"         );
    iEvent.put(std::move( vector_mus_gfit_nlayers3D       	      ), branchprefix_ + "gfitnlayers3D"       );
    iEvent.put(std::move( vector_mus_gfit_nlayersLost     	      ), branchprefix_ + "gfitnlayersLost"     );
    iEvent.put(std::move( vector_mus_gfit_validPixelHits  	      ), branchprefix_ + "gfitvalidPixelHits"  );
    iEvent.put(std::move( vector_mus_gfit_lostPixelHits   	      ), branchprefix_ + "gfitlostPixelHits"   );
    iEvent.put(std::move( vector_mus_gfit_exp_innerlayers 	      ), branchprefix_ + "gfitexpinnerlayer"  );
    iEvent.put(std::move( vector_mus_gfit_exp_outerlayers 	      ), branchprefix_ + "gfitexpouterlayer"  );
    iEvent.put(std::move( vector_mus_gfit_lostHits        	      ), branchprefix_ + "gfitlostHits"        );
    iEvent.put(std::move( vector_mus_gfit_etaErr          	      ), branchprefix_ + "gfitetaErr"          );
    iEvent.put(std::move( vector_mus_gfit_phiErr          	      ), branchprefix_ + "gfitphiErr"          );
    iEvent.put(std::move( vector_mus_gfit_trk_charge              ), branchprefix_ + "gfittrkcharge"       );

    ////////////                                                                       
    // Best   //
    ////////////

    iEvent.put(std::move( vector_mus_bfit_d0                      ), branchprefix_ + "bfitd0"             );
    iEvent.put(std::move( vector_mus_bfit_z0                      ), branchprefix_ + "bfitz0"             );
    iEvent.put(std::move( vector_mus_bfit_d0Err                   ), branchprefix_ + "bfitd0Err"          );
    iEvent.put(std::move( vector_mus_bfit_z0Err                   ), branchprefix_ + "bfitz0Err"          );
    iEvent.put(std::move( vector_mus_bfit_d0corr                  ), branchprefix_ + "bfitd0corr"         );
    iEvent.put(std::move( vector_mus_bfit_z0corr                  ), branchprefix_ + "bfitz0corr"         );
    iEvent.put(std::move( vector_mus_bfit_qoverp                  ), branchprefix_ + "bfitqoverp"         );
    iEvent.put(std::move( vector_mus_bfit_qoverpError             ), branchprefix_ + "bfitqoverpError"    );
    iEvent.put(std::move( vector_mus_bfit_chi2                    ), branchprefix_ + "bfitchi2"           );
    iEvent.put(std::move( vector_mus_bfit_ndof                    ), branchprefix_ + "bfitndof"           );
    iEvent.put(std::move( vector_mus_bfit_validHits               ), branchprefix_ + "bfitvalidHits"      );
    iEvent.put(std::move( vector_mus_bfit_validSTAHits            ), branchprefix_ + "bfitvalidSTAHits"   );
    iEvent.put(std::move( vector_mus_bfit_validSiHits             ), branchprefix_ + "bfitvalidSiHits"    );

    iEvent.put(std::move( vector_mus_bfit_vertex_p4               ), branchprefix_ + "bfitvertexp4"             );
    iEvent.put(std::move( vector_mus_bfit_d0corrPhi               ), branchprefix_ + "bfitd0corrPhi"       );
    iEvent.put(std::move( vector_mus_bfit_d0phiCov        	      ), branchprefix_ + "bfitd0phiCov"        );
    iEvent.put(std::move( vector_mus_bfit_qualityMask     	      ), branchprefix_ + "bfitqualityMask"     );
    iEvent.put(std::move( vector_mus_bfit_nlayers         	      ), branchprefix_ + "bfitnlayers"         );
    iEvent.put(std::move( vector_mus_bfit_nlayers3D       	      ), branchprefix_ + "bfitnlayers3D"       );
    iEvent.put(std::move( vector_mus_bfit_nlayersLost     	      ), branchprefix_ + "bfitnlayersLost"     );
    iEvent.put(std::move( vector_mus_bfit_validPixelHits  	      ), branchprefix_ + "bfitvalidPixelHits"  );
    iEvent.put(std::move( vector_mus_bfit_lostPixelHits   	      ), branchprefix_ + "bfitlostPixelHits"   );
    iEvent.put(std::move( vector_mus_bfit_exp_innerlayers 	      ), branchprefix_ + "bfitexpinnerlayer"  );
    iEvent.put(std::move( vector_mus_bfit_exp_outerlayers 	      ), branchprefix_ + "bfitexpouterlayer"  );
    iEvent.put(std::move( vector_mus_bfit_lostHits        	      ), branchprefix_ + "bfitlostHits"        );
    iEvent.put(std::move( vector_mus_bfit_etaErr          	      ), branchprefix_ + "bfitetaErr"          );
    iEvent.put(std::move( vector_mus_bfit_phiErr          	      ), branchprefix_ + "bfitphiErr"          );
    iEvent.put(std::move( vector_mus_bfit_trk_charge              ), branchprefix_ + "bfittrkcharge"       );

    //////////////////
    // Muon Quality //
    //////////////////

    iEvent.put(std::move( vector_mus_updatedSta         ), branchprefix_ + "updatedSta"         );
    iEvent.put(std::move( vector_mus_tightMatch         ), branchprefix_ + "tightMatch"         );
    iEvent.put(std::move( vector_mus_glbKink            ), branchprefix_ + "glbKink"            );
    iEvent.put(std::move( vector_mus_trkRelChi2         ), branchprefix_ + "trkRelChi2"         );
    iEvent.put(std::move( vector_mus_staRelChi2         ), branchprefix_ + "staRelChi2"         );
    iEvent.put(std::move( vector_mus_localDistance      ), branchprefix_ + "localDistance"      );
    iEvent.put(std::move( vector_mus_globalDeltaEtaPhi  ), branchprefix_ + "globalDeltaEtaPhi"  );
    iEvent.put(std::move( vector_mus_glbTrackProbability), branchprefix_ + "glbTrackProbability");


    ////////////
    // Timing //
    ////////////

    iEvent.put(std::move( vector_mus_timeNumStationsUsed          ), branchprefix_ + "timeNumStationsUsed"      );
    iEvent.put(std::move( vector_mus_timeAtIpInOut                ), branchprefix_ + "timeAtIpInOut"            );
    iEvent.put(std::move( vector_mus_timeAtIpInOutErr             ), branchprefix_ + "timeAtIpInOutErr"         );
    iEvent.put(std::move( vector_mus_timeAtIpOutIn                ), branchprefix_ + "timeAtIpOutIn"            );
    iEvent.put(std::move( vector_mus_timeAtIpOutInErr             ), branchprefix_ + "timeAtIpOutInErr"         );
    iEvent.put(std::move( vector_mus_timeDirection                ), branchprefix_ + "timeDirection"            );


    ///////////
    // Muons //
    ///////////

    iEvent.put(std::move( vector_mus_goodmask                ), branchprefix_ + "goodmask"                );
    iEvent.put(std::move( vector_mus_nmatches                ), branchprefix_ + "nmatches"                );
    iEvent.put(std::move( vector_mus_nOverlaps               ), branchprefix_ + "nOverlaps"               );
    iEvent.put(std::move( vector_mus_overlap0                ), branchprefix_ + "overlap0"                );
    iEvent.put(std::move( vector_mus_overlap1                ), branchprefix_ + "overlap1"                );
    iEvent.put(std::move( vector_mus_mass                    ), branchprefix_ + "mass"                    );
    iEvent.put(std::move( vector_mus_isRPCMuon               ), branchprefix_ + "isRPCMuon"               );
    

    ////////////
    // Energy //
    ////////////

    iEvent.put(std::move( vector_mus_e_em               ), branchprefix_ + "eem"                );
    iEvent.put(std::move( vector_mus_e_had              ), branchprefix_ + "ehad"               );
    iEvent.put(std::move( vector_mus_e_ho               ), branchprefix_ + "eho"                );
    iEvent.put(std::move( vector_mus_e_emS9             ), branchprefix_ + "eemS9"              );
    iEvent.put(std::move( vector_mus_e_hadS9            ), branchprefix_ + "ehadS9"             );
    iEvent.put(std::move( vector_mus_e_hoS9             ), branchprefix_ + "ehoS9"              );
    iEvent.put(std::move( vector_mus_emS25              ), branchprefix_ + "emS25"              );
    iEvent.put(std::move( vector_mus_emMax              ), branchprefix_ + "emMax"              );
    iEvent.put(std::move( vector_mus_hadMax             ), branchprefix_ + "hadMax"             );
    iEvent.put(std::move( vector_mus_ecal_rawId         ), branchprefix_ + "ecalrawId"          );
    iEvent.put(std::move( vector_mus_hcal_rawId         ), branchprefix_ + "hcalrawId"          );
  
    ///////////////
    // Isolation //
    ///////////////

    iEvent.put(std::move( vector_mus_iso03_hoEt         ), branchprefix_ + "iso03hoEt"          );
    iEvent.put(std::move( vector_mus_iso03_njets        ), branchprefix_ + "iso03njets"         );
    iEvent.put(std::move( vector_mus_iso05_sumPt        ), branchprefix_ + "iso05sumPt"         );
    iEvent.put(std::move( vector_mus_iso05_emEt         ), branchprefix_ + "iso05emEt"          );
    iEvent.put(std::move( vector_mus_iso05_hadEt        ), branchprefix_ + "iso05hadEt"         );
    iEvent.put(std::move( vector_mus_iso05_hoEt         ), branchprefix_ + "iso05hoEt"          );
    iEvent.put(std::move( vector_mus_iso05_ntrk         ), branchprefix_ + "iso05ntrk"          );


    ////////////
    // Tracks //
    ////////////

    iEvent.put(std::move( vector_mus_muonBestTrackType  ), branchprefix_ + "muonBestTrackType"  );
    iEvent.put(std::move( vector_mus_vertex_p4          ), branchprefix_ + "vertexp4"           );
    iEvent.put(std::move( vector_mus_d0                 ), branchprefix_ + "d0"                 );
    iEvent.put(std::move( vector_mus_z0                 ), branchprefix_ + "z0"                 );
    iEvent.put(std::move( vector_mus_d0corr             ), branchprefix_ + "d0corr"             );
    iEvent.put(std::move( vector_mus_z0corr             ), branchprefix_ + "z0corr"             );
    iEvent.put(std::move( vector_mus_vertexphi          ), branchprefix_ + "vertexphi"          );
    iEvent.put(std::move( vector_mus_chi2               ), branchprefix_ + "chi2"               );
    iEvent.put(std::move( vector_mus_ndof               ), branchprefix_ + "ndof"               );
    iEvent.put(std::move( vector_mus_etaErr             ), branchprefix_ + "etaErr"             );
    iEvent.put(std::move( vector_mus_phiErr             ), branchprefix_ + "phiErr"             );
    iEvent.put(std::move( vector_mus_trk_charge         ), branchprefix_ + "trkcharge"          );
    iEvent.put(std::move( vector_mus_qoverp             ), branchprefix_ + "qoverp"             );
    iEvent.put(std::move( vector_mus_qoverpError        ), branchprefix_ + "qoverpError"        );
    iEvent.put(std::move( vector_mus_d0corrPhi          ), branchprefix_ + "d0corrPhi"      	   );
    iEvent.put(std::move( vector_mus_d0phiCov           ), branchprefix_ + "d0phiCov"       	   );
    iEvent.put(std::move( vector_mus_qualityMask        ), branchprefix_ + "qualityMask"    	   );
    iEvent.put(std::move( vector_mus_nlayers3D          ), branchprefix_ + "nlayers3D"      	   );
    iEvent.put(std::move( vector_mus_nlayersLost        ), branchprefix_ + "nlayersLost"    	   );
    iEvent.put(std::move( vector_mus_lostPixelHits      ), branchprefix_ + "lostPixelHits"  	   );
  
    ////////                                  
    // PF //
    ////////

    iEvent.put(std::move( vector_mus_isoMeanDRR03_pf_ChargedHadronPt                ), branchprefix_ + "isoMeanDRR03pfChargedHadronPt"             );
    iEvent.put(std::move( vector_mus_isoMeanDRR03_pf_ChargedParticlePt              ), branchprefix_ + "isoMeanDRR03pfChargedParticlePt"           );
    iEvent.put(std::move( vector_mus_isoMeanDRR03_pf_NeutralHadronEt                ), branchprefix_ + "isoMeanDRR03pfNeutralHadronEt"             );
    iEvent.put(std::move( vector_mus_isoMeanDRR03_pf_PhotonEt                       ), branchprefix_ + "isoMeanDRR03pfPhotonEt"                    );
    iEvent.put(std::move( vector_mus_isoMeanDRR03_pf_sumNeutralHadronEtHighThreshold), branchprefix_ + "isoMeanDRR03pfNeutralHadronEtHighThreshold");
    iEvent.put(std::move( vector_mus_isoMeanDRR03_pf_sumPhotonEtHighThreshold       ), branchprefix_ + "isoMeanDRR03pfPhotonEtHighThreshold"       );
    iEvent.put(std::move( vector_mus_isoMeanDRR03_pf_PUPt                           ), branchprefix_ + "isoMeanDRR03pfPUPt"                        );

    iEvent.put(std::move( vector_mus_isoMeanDRR04_pf_ChargedHadronPt                ), branchprefix_ + "isoMeanDRR04pfChargedHadronPt"             );
    iEvent.put(std::move( vector_mus_isoMeanDRR04_pf_ChargedParticlePt              ), branchprefix_ + "isoMeanDRR04pfChargedParticlePt"           );
    iEvent.put(std::move( vector_mus_isoMeanDRR04_pf_NeutralHadronEt                ), branchprefix_ + "isoMeanDRR04pfNeutralHadronEt"             );
    iEvent.put(std::move( vector_mus_isoMeanDRR04_pf_PhotonEt                       ), branchprefix_ + "isoMeanDRR04pfPhotonEt"                    );
    iEvent.put(std::move( vector_mus_isoMeanDRR04_pf_sumNeutralHadronEtHighThreshold), branchprefix_ + "isoMeanDRR04pfNeutralHadronEtHighThreshold");
    iEvent.put(std::move( vector_mus_isoMeanDRR04_pf_sumPhotonEtHighThreshold       ), branchprefix_ + "isoMeanDRR04pfPhotonEtHighThreshold"       );
    iEvent.put(std::move( vector_mus_isoMeanDRR04_pf_PUPt                           ), branchprefix_ + "isoMeanDRR04pfPUPt"                        );

    iEvent.put(std::move( vector_mus_isoSumDRR03_pf_ChargedHadronPt                ), branchprefix_ + "isoSumDRR03pfChargedHadronPt"             );
    iEvent.put(std::move( vector_mus_isoSumDRR03_pf_ChargedParticlePt              ), branchprefix_ + "isoSumDRR03pfChargedParticlePt"           );
    iEvent.put(std::move( vector_mus_isoSumDRR03_pf_NeutralHadronEt                ), branchprefix_ + "isoSumDRR03pfNeutralHadronEt"             );
    iEvent.put(std::move( vector_mus_isoSumDRR03_pf_PhotonEt                       ), branchprefix_ + "isoSumDRR03pfPhotonEt"                    );
    iEvent.put(std::move( vector_mus_isoSumDRR03_pf_sumNeutralHadronEtHighThreshold), branchprefix_ + "isoSumDRR03pfNeutralHadronEtHighThreshold");
    iEvent.put(std::move( vector_mus_isoSumDRR03_pf_sumPhotonEtHighThreshold       ), branchprefix_ + "isoSumDRR03pfPhotonEtHighThreshold"       );
    iEvent.put(std::move( vector_mus_isoSumDRR03_pf_PUPt                           ), branchprefix_ + "isoSumDRR03pfPUPt"                        );

    iEvent.put(std::move( vector_mus_isoSumDRR04_pf_ChargedHadronPt                ), branchprefix_ + "isoSumDRR04pfChargedHadronPt"             );
    iEvent.put(std::move( vector_mus_isoSumDRR04_pf_ChargedParticlePt              ), branchprefix_ + "isoSumDRR04pfChargedParticlePt"           );
    iEvent.put(std::move( vector_mus_isoSumDRR04_pf_NeutralHadronEt                ), branchprefix_ + "isoSumDRR04pfNeutralHadronEt"             );
    iEvent.put(std::move( vector_mus_isoSumDRR04_pf_PhotonEt                       ), branchprefix_ + "isoSumDRR04pfPhotonEt"                    );
    iEvent.put(std::move( vector_mus_isoSumDRR04_pf_sumNeutralHadronEtHighThreshold), branchprefix_ + "isoSumDRR04pfNeutralHadronEtHighThreshold");
    iEvent.put(std::move( vector_mus_isoSumDRR04_pf_sumPhotonEtHighThreshold       ), branchprefix_ + "isoSumDRR04pfPhotonEtHighThreshold"       );
    iEvent.put(std::move( vector_mus_isoSumDRR04_pf_PUPt                           ), branchprefix_ + "isoSumDRR04pfPUPt"                        );

    ///////////
    // IP 3D //
    ///////////

    iEvent.put(std::move( vector_mus_bs3d                         ), branchprefix_ + "bs3d"               );
    iEvent.put(std::move( vector_mus_bs3derr                      ), branchprefix_ + "bs3derr"            );
    iEvent.put(std::move( vector_mus_bs2d                         ), branchprefix_ + "bs2d"               );
    iEvent.put(std::move( vector_mus_bs2derr                      ), branchprefix_ + "bs2derr"            );

} //

//define this as a plug-in
DEFINE_FWK_MODULE(MuonExtraMaker);
