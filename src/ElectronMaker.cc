//-*- C++ -*- 
//
// Package:    ElectronMaker
// Class:      ElectronMaker
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronMaker.cc,v 1.89 2012/08/16 00:00:27 slava77 Exp $


//System include files
#include <memory>
#include <math.h>

//User include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CMS3/NtupleMaker/interface/ElectronMaker.h"
#include "CMS3/NtupleMaker/interface/MatchUtilities.h"
#include "CMS3/NtupleMaker/interface/MCUtilities.h"
#include "CMS3/NtupleMaker/interface/EgammaFiduciality.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Math/VectorUtil.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"

using namespace reco;
using namespace edm;
using namespace std;

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
typedef Ref<edmNew::DetSetVector<SiStripCluster>,SiStripCluster > ClusterRef;
typedef Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster > pixel_ClusterRef;

// constructors and destructor
ElectronMaker::ElectronMaker(const ParameterSet& iConfig) {

    //get setup parameters
    electronVetoIdMapToken_   = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"));
    electronLooseIdMapToken_  = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"));
    electronMediumIdMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"));
    electronTightIdMapToken_  = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"));

    electronsInputTag_        = iConfig.getParameter<edm::InputTag> ("electronsInputTag"        );
    beamSpotInputTag_         = iConfig.getParameter<edm::InputTag> ("beamSpotInputTag"         );
    trksInputTag_             = iConfig.getParameter<edm::InputTag> ("trksInputTag"             );
    gsftracksInputTag_        = iConfig.getParameter<edm::InputTag> ("gsftracksInputTag"        );
    cms2scsseeddetidInputTag_ = iConfig.getParameter<edm::InputTag> ("cms2scsseeddetidInputTag" );
    eidLHTag_                 = iConfig.getParameter<edm::InputTag> ("eidLHTag"                 );
    pfCandsInputTag           = iConfig.getParameter<edm::InputTag> ("pfCandsInputTag"          );
    vtxInputTag               = iConfig.getParameter<edm::InputTag> ("vtxInputTag"              );
    // pfIsoCharged03InputTag    = iConfig.getParameter<edm::InputTag> ("pfIsoCharged03InputTag"   );
    // pfIsoGamma03InputTag      = iConfig.getParameter<edm::InputTag> ("pfIsoGamma03InputTag"     );
    // pfIsoNeutral03InputTag    = iConfig.getParameter<edm::InputTag> ("pfIsoNeutral03InputTag"   );
    // pfIsoCharged04InputTag    = iConfig.getParameter<edm::InputTag> ("pfIsoCharged04InputTag"   );
    // pfIsoGamma04InputTag      = iConfig.getParameter<edm::InputTag> ("pfIsoGamma04InputTag"     );
    // pfIsoNeutral04InputTag    = iConfig.getParameter<edm::InputTag> ("pfIsoNeutral04InputTag"   );
    

    recoConversionInputTag_   = iConfig.getParameter<edm::InputTag> ("recoConversionInputTag"   );
    rhoInputTag_              = iConfig.getParameter<edm::InputTag> ("rhoInputTag"              );
    beamSpot_tag_             = iConfig.getParameter<edm::InputTag> ("beamSpotTag"              );

    minAbsDist_               = iConfig.getParameter<double>          ("minAbsDist"              );
    minAbsDcot_               = iConfig.getParameter<double>          ("minAbsDcot"              );
    minSharedFractionOfHits_  = iConfig.getParameter<double>          ("minSharedFractionOfHits" );
    aliasprefix_              = iConfig.getUntrackedParameter<string> ("aliasPrefix"             );

    mtsTransform_ = 0;
    clusterTools_ = 0;
 

    produces<unsigned int>       ("evtnels"                    ).setBranchAlias("evt_nels"                   ); //number of electrons in event

    // ECAL related (superCluster) variables
    produces<vector<int> >       ("elsnSeed"                   ).setBranchAlias("els_nSeed"                  );
    produces<vector<float> >     ("elseSC"                     ).setBranchAlias("els_eSC"                    );
    produces<vector<float> >     ("elsetaSC"                   ).setBranchAlias("els_etaSC"                  );
    produces<vector<float> >     ("elsphiSC"                   ).setBranchAlias("els_phiSC"                  );
    produces<vector<float> >     ("elseSCRaw"                  ).setBranchAlias("els_eSCRaw"                 );
    produces<vector<float> >     ("elseSCPresh"                ).setBranchAlias("els_eSCPresh"               );
    produces<vector<int> >       ("elsfiduciality"             ).setBranchAlias("els_fiduciality"            );
    produces<vector<int> >       ("elstype"                    ).setBranchAlias("els_type"                   );
//    produces<vector<int> >       ("elsscindex"                 ).setBranchAlias("els_scindex"                );
    produces<vector<float> >     ("elsetaSCwidth"              ).setBranchAlias("els_etaSCwidth"             );
    produces<vector<float> >     ("elsphiSCwidth"              ).setBranchAlias("els_phiSCwidth"             );

    // Corrections and uncertainties
    //
    produces<vector<float> >     ("elsecalEnergy"              ).setBranchAlias("els_ecalEnergy"             );
    produces<vector<float> >     ("elsecalEnergyError"         ).setBranchAlias("els_ecalEnergyError"        );
    produces<vector<float> >     ("elstrackMomentumError"      ).setBranchAlias("els_trackMomentumError"     );

    // ID variables
    //
    //produces<vector<float> >     ("elslh"                      ).setBranchAlias("els_lh"                     );
    produces<vector<float> >     ("elsmva"                     ).setBranchAlias("els_mva"                    );

    produces<vector<float> >     ("elsdEtaIn"                  ).setBranchAlias("els_dEtaIn"                 );
    produces<vector<float> >     ("elsdEtaOut"                 ).setBranchAlias("els_dEtaOut"                );
    produces<vector<float> >     ("elsdPhiIn"                  ).setBranchAlias("els_dPhiIn"                 );
    produces<vector<float> >     ("elsdPhiOut"                 ).setBranchAlias("els_dPhiOut"                );
    produces<vector<float> >     ("elsdPhiInPhiOut"            ).setBranchAlias("els_dPhiInPhiOut"           );
    produces<vector<float> >     ("elsfbrem"                   ).setBranchAlias("els_fbrem"                  );
    produces<vector<float> >     ("elseSeed"                   ).setBranchAlias("els_eSeed"                  );
    produces<vector<float> >     ("elseOverPIn"                ).setBranchAlias("els_eOverPIn"               );
    produces<vector<float> >     ("elseSeedOverPOut"           ).setBranchAlias("els_eSeedOverPOut"          );
    produces<vector<float> >     ("elseSeedOverPIn"            ).setBranchAlias("els_eSeedOverPIn"           );
    produces<vector<float> >     ("elseOverPOut"               ).setBranchAlias("els_eOverPOut"              );

    produces<vector<float> >     ("elshOverE"                  ).setBranchAlias("els_hOverE"                 );
    produces<vector<float> >     ("elshOverEBC"                ).setBranchAlias("els_hOverEBC"               );
    produces<vector<float> >     ("elshcalDepth1OverEcal"      ).setBranchAlias("els_hcalDepth1OverEcal"     );
    produces<vector<float> >     ("elshcalDepth2OverEcal"      ).setBranchAlias("els_hcalDepth2OverEcal"     );

    //produces<vector<float> >     ("elssigmaPhiPhi"             ).setBranchAlias("els_sigmaPhiPhi"            );
    produces<vector<float> >     ("elssigmaIPhiIPhi"           ).setBranchAlias("els_sigmaIPhiIPhi"          );
    //produces<vector<float> >     ("elssigmaIEtaIPhi"           ).setBranchAlias("els_sigmaIEtaIPhi"          );
    produces<vector<float> >     ("elssigmaEtaEta"             ).setBranchAlias("els_sigmaEtaEta"            );
    produces<vector<float> >     ("elssigmaIEtaIEta"           ).setBranchAlias("els_sigmaIEtaIEta"          );
    //produces<vector<float> >     ("elssigmaIPhiIPhiSC"         ).setBranchAlias("els_sigmaIPhiIPhiSC"        );
    //produces<vector<float> >     ("elssigmaIEtaIEtaSC"         ).setBranchAlias("els_sigmaIEtaIEtaSC"        );

    produces<vector<float> >     ("else2x5Max"                 ).setBranchAlias("els_e2x5Max"                );
    produces<vector<float> >     ("else1x5"                    ).setBranchAlias("els_e1x5"                   );
    produces<vector<float> >     ("else5x5"                    ).setBranchAlias("els_e5x5"                   );
//    produces<vector<float> >     ("else3x3"                    ).setBranchAlias("els_e3x3"                   );
//    produces<vector<float> >     ("elseMax"                    ).setBranchAlias("els_eMax"                   );

    produces<vector<float> >     ("elsdeltaEtaEleClusterTrackAtCalo").setBranchAlias("els_deltaEtaEleClusterTrackAtCalo");
    produces<vector<float> >     ("elsdeltaPhiEleClusterTrackAtCalo").setBranchAlias("els_deltaPhiEleClusterTrackAtCalo");
    produces<vector<bool > >     ("elsisGsfCtfScPixChargeConsistent").setBranchAlias("els_isGsfCtfScPixChargeConsistent");

    // predefined ID decisions
    // http://cmslxr.fnal.gov/lxr/source/DataFormats/EgammaCandidates/interface/GsfElectron.h
    produces<vector<int> >       ("elsclass"                   ).setBranchAlias("els_class"                  );

    // Phys 14 predefined ID decisions
    produces<vector<int> >       ("passVetoId"                 ).setBranchAlias("els_passVetoId"                 );
    produces<vector<int> >       ("passLooseId"                ).setBranchAlias("els_passLooseId"                );
    produces<vector<int> >       ("passMediumId"               ).setBranchAlias("els_passMediumId"               );
    produces<vector<int> >       ("passTightId"                ).setBranchAlias("els_passTightId"                );

    // for the ID definitions, see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideElectronID
    // the decisions should be the SAME as the els_pat_*id branches made by PATElectronMaker
    produces<vector<int> >       ("elscategory"                ).setBranchAlias("els_category"               );

    // isolation variables
    //
    produces<vector<float> >     ("elstkIso"                  ).setBranchAlias("els_tkIso"                  );
    produces<vector<float> >     ("elsecalIso"                ).setBranchAlias("els_ecalIso"                );
    produces<vector<float> >     ("elshcalIso"                ).setBranchAlias("els_hcalIso"                );
    produces<vector<float> >     ("elshcalDepth1TowerSumEt"   ).setBranchAlias("els_hcalDepth1TowerSumEt"   );
    produces<vector<float> >     ("elshcalDepth2TowerSumEt"   ).setBranchAlias("els_hcalDepth2TowerSumEt"   );

    produces<vector<float> >     ("elstkIso04"                ).setBranchAlias("els_tkIso04"                );
    produces<vector<float> >     ("elsecalIso04"              ).setBranchAlias("els_ecalIso04"              );
    produces<vector<float> >     ("elshcalIso04"              ).setBranchAlias("els_hcalIso04"              );
    produces<vector<float> >     ("elshcalDepth1TowerSumEt04" ).setBranchAlias("els_hcalDepth1TowerSumEt04" );
    produces<vector<float> >     ("elshcalDepth2TowerSumEt04" ).setBranchAlias("els_hcalDepth2TowerSumEt04" );
//    produces<vector<float> >     ("elsiso03pf"                ).setBranchAlias("els_iso03_pf"               ); // pf isolation in cone of 0.3, 1 GeV threshold
//    produces<vector<float> >     ("elsiso04pf"                ).setBranchAlias("els_iso04_pf"               ); // pf isolation in cone of 0.4, 1 GeV threshold

//    produces<vector<float> >     ("elsiso03pfch"              ).setBranchAlias("els_iso03_pf_ch"            ); // pf isolation in cone of 0.3, charged only
//    produces<vector<float> >     ("elsiso03pfgamma05"         ).setBranchAlias("els_iso03_pf_gamma05"       ); // pf isolation in cone of 0.3, photons only with threshold 0.5 GeV
//    produces<vector<float> >     ("elsiso03pfnhad05"          ).setBranchAlias("els_iso03_pf_nhad05"        ); // pf isolation in cone of 0.3, neutral hadrons only with threshold 0.5 GeV
 //   produces<vector<float> >     ("elsiso04pfch"              ).setBranchAlias("els_iso04_pf_ch"            ); // pf isolation in cone of 0.3, charged only
//    produces<vector<float> >     ("elsiso04pfgamma05"         ).setBranchAlias("els_iso04_pf_gamma05"       ); // pf isolation in cone of 0.3, photons only with threshold 0.5 GeV
//    produces<vector<float> >     ("elsiso04pfnhad05"          ).setBranchAlias("els_iso04_pf_nhad05"        ); // pf isolation in cone of 0.3, neutral hadrons only with threshold 0.5 GeV


    // pf isolation variables
    produces<vector<float> >     ("elspfChargedHadronIso").setBranchAlias("els_pfChargedHadronIso");
    produces<vector<float> >     ("elspfNeutralHadronIso").setBranchAlias("els_pfNeutralHadronIso");
    produces<vector<float> >     ("elspfPhotonIso"       ).setBranchAlias("els_pfPhotonIso"       );
    produces<vector<float> >     ("elspfPUIso"           ).setBranchAlias("els_pfPUIso"           );


    // track variables
    //
    produces<vector<int> >       ("elscharge"        ).setBranchAlias("els_charge"         ); //candidate charge
    produces<vector<int> >       ("elssccharge"      ).setBranchAlias("els_sccharge"       );
    produces<vector<int> >       ("elstrkcharge"     ).setBranchAlias("els_trk_charge"     );
    produces<vector<float> >     ("elschi2"          ).setBranchAlias("els_chi2"           );
    produces<vector<float> >     ("elsndof"          ).setBranchAlias("els_ndof"           );
    produces<vector<float> >     ("elsckfchi2"       ).setBranchAlias("els_ckf_chi2"       );
    produces<vector<float> >     ("elsckfndof"       ).setBranchAlias("els_ckf_ndof"       );
    produces<vector<int> >       ("elsvalidHits"     ).setBranchAlias("els_validHits"      ); //number of used hits in fit
    produces<vector<int> >       ("elslostHits"      ).setBranchAlias("els_lostHits"       ); //number of lost hits in fit
    produces<vector<float> >     ("elsd0"            ).setBranchAlias("els_d0"             );
    produces<vector<float> >     ("elsz0"            ).setBranchAlias("els_z0"             );
    produces<vector<float> >     ("elsdxyPV"         ).setBranchAlias("els_dxyPV"          );
    produces<vector<float> >     ("elsdzPV"          ).setBranchAlias("els_dzPV"           );
    produces<vector<float> >     ("elsd0Err"         ).setBranchAlias("els_d0Err"          );
    produces<vector<float> >     ("elsz0Err"         ).setBranchAlias("els_z0Err"          );
    produces<vector<float> >     ("elsd0corr"        ).setBranchAlias("els_d0corr"         );
    produces<vector<float> >     ("elsd0corrPhi"     ).setBranchAlias("els_d0corrPhi"      );
    produces<vector<float> >     ("elsd0phiCov"      ).setBranchAlias("els_d0phiCov"       );
    produces<vector<float> >     ("elsz0corr"        ).setBranchAlias("els_z0corr"         );
    produces<vector<float> >     ("elsptErr"         ).setBranchAlias("els_ptErr"          );
    produces<vector<float> >     ("elsptErrGsf"      ).setBranchAlias("els_ptErrGsf"       );
    produces<vector<float> >     ("elsetaErr"        ).setBranchAlias("els_etaErr"         );
    produces<vector<float> >     ("elsphiErr"        ).setBranchAlias("els_phiErr"         );
    produces<vector<float> >     ("elsip3d"          ).setBranchAlias("els_ip3d"           ); // Ip3d from normal vertex
    produces<vector<float> >     ("elsip3derr"       ).setBranchAlias("els_ip3derr"        ); // Ip3d error from normal vertex
    produces<vector<float> >     ("elsip2d"          ).setBranchAlias("els_ip2d"           ); // Ip2d from normal vertex
    produces<vector<float> >     ("elsip2derr"       ).setBranchAlias("els_ip2derr"        ); // Ip2d error from normal vertex
    produces<vector<float> >     ("elsbs3d"          ).setBranchAlias("els_bs3d"           ); // Ip3d from normal vertex
    produces<vector<float> >     ("elsbs3derr"       ).setBranchAlias("els_bs3derr"        ); // Ip3d error from normal vertex
    produces<vector<float> >     ("elsbs2d"          ).setBranchAlias("els_bs2d"           ); // Ip2d from normal vertex
    produces<vector<float> >     ("elsbs2derr"       ).setBranchAlias("els_bs2derr"        ); // Ip2d error from normal vertex
    produces<vector<int> >       ("elsckflaywithmeas").setBranchAlias("els_ckf_laywithmeas");
    produces<vector<int> >       ("elsckfcharge"     ).setBranchAlias("els_ckf_charge"     );

    // LorentzVectors
    //
    produces<vector<LorentzVector> >  ("elsp4"    ).setBranchAlias("els_p4"     );
    produces<vector<LorentzVector> >  ("elstrkp4" ).setBranchAlias("els_trk_p4" );
    produces<vector<LorentzVector> >  ("elsp4In"  ).setBranchAlias("els_p4In"   );
    produces<vector<LorentzVector> >  ("elsp4Out" ).setBranchAlias("els_p4Out"  );

    // Vertex
    //
    produces<vector<LorentzVector> >  ("elsvertexp4").setBranchAlias("els_vertex_p4");
    produces<vector<LorentzVector> >  ("elstrkvertexp4").setBranchAlias("els_trk_vertex_p4");

    //Hit Pattern information
    //
    produces<vector<LorentzVector> >  ("elsinnerposition"  ).setBranchAlias("els_inner_position"  );
    produces<vector<LorentzVector> >  ("elsouterposition"  ).setBranchAlias("els_outer_position"  );
    produces<vector<int> >            ("elsvalidpixelhits" ).setBranchAlias("els_valid_pixelhits" );
    produces<vector<int> >            ("elslostpixelhits"  ).setBranchAlias("els_lost_pixelhits"  );
    produces<vector<int> >            ("elsnlayers"        ).setBranchAlias("els_nlayers"         );
    produces<vector<int> >            ("elsnlayers3D"      ).setBranchAlias("els_nlayers3D"       );
    produces<vector<int> >            ("elsnlayersLost"    ).setBranchAlias("els_nlayersLost"     );
    //produces<vector<int> >            ("elslayer1sizerphi" ).setBranchAlias("els_layer1_sizerphi" ); 
    //produces<vector<int> >            ("elslayer1sizerz"   ).setBranchAlias("els_layer1_sizerz"   ); 
    //produces<vector<float> >          ("elslayer1charge"   ).setBranchAlias("els_layer1_charge"   ); 
    //produces<vector<int> >            ("elslayer1det"      ).setBranchAlias("els_layer1_det"      );
    //produces<vector<int> >            ("elslayer1layer"    ).setBranchAlias("els_layer1_layer"    ); 
    produces<vector<int> >            ("elsexpinnerlayers" ).setBranchAlias("els_exp_innerlayers" );
    produces<vector<int> >            ("elsexpouterlayers" ).setBranchAlias("els_exp_outerlayers" );   

    //CTF track matching stuff
    produces<vector<float>  >    ("elstrkshFrac" ).setBranchAlias("els_trkshFrac" );
    produces<vector<float>  >    ("elstrkdr"     ).setBranchAlias("els_trkdr"     );


    //These are vectors of vectors, holding the candidate conversion partners. 
    //produces<vector<vector<float> >   >       ("elsconvsdist"        ).setBranchAlias("els_convs_dist"        );
    //produces<vector<vector<float> >   >       ("elsconvsdcot"        ).setBranchAlias("els_convs_dcot"        );
    //produces<vector<vector<float> >   >       ("elsconvsradius"      ).setBranchAlias("els_convs_radius"      );  //signed radius of conversion
    //produces<vector<vector<int>   >   >       ("elsconvstkidx"       ).setBranchAlias("els_convs_tkidx"       );  //index of partner track
    //produces<vector<vector<int>   >   >       ("elsconvsgsftkidx"    ).setBranchAlias("els_convs_gsftkidx"    );  //index of the GSF partner track, if thats where we find it
    //produces<vector<vector<int>   >   >       ("elsconvsdelMissHits" ).setBranchAlias("els_convs_delMissHits" );  //Delta Missing Hits between the electron and partner track
    //produces<vector<vector<int>   >   >       ("elsconvsflag"        ).setBranchAlias("els_convs_flag"        );


    //conversion stuff - This is the "best" conversion partner, defined as the one that 
    //has the minimum sqrt(dist*dist + dcot*dcot)
    //produces<vector<LorentzVector> >  ("elsconvposp4"       ).setBranchAlias("els_conv_pos_p4"      );  //position of conversion
    //produces<vector<float>    >       ("elsconvdist"        ).setBranchAlias("els_conv_dist"        );
    //produces<vector<float>    >       ("elsconvdcot"        ).setBranchAlias("els_conv_dcot"        );
    //produces<vector<float>    >       ("elsconvradius"      ).setBranchAlias("els_conv_radius"      );  //signed radius of conversion
    //produces<vector<int>      >       ("elsconvtkidx"       ).setBranchAlias("els_conv_tkidx"       );  //index of partner track
    //produces<vector<int>      >       ("elsconvgsftkidx"    ).setBranchAlias("els_conv_gsftkidx"    );  //index of the GSF partner track, if thats where we find it
    //produces<vector<int>      >       ("elsconvdelMissHits" ).setBranchAlias("els_conv_delMissHits" );  //Delta Missing Hits between the electron and partner track
    //produces<vector<int>      >       ("elsconvflag"        ).setBranchAlias("els_conv_flag"        );
  
    //produces<vector<float>    >       ("elsconvolddist"        ).setBranchAlias("els_conv_old_dist"        );
    //produces<vector<float>    >       ("elsconvolddcot"        ).setBranchAlias("els_conv_old_dcot"        );
    //produces<vector<float>    >       ("elsconvoldradius"      ).setBranchAlias("els_conv_old_radius"      );  //signed radius of conversion
    //produces<vector<int>      >       ("elsconvoldtkidx"       ).setBranchAlias("els_conv_old_tkidx"       );  //index of partner track
    //produces<vector<int>      >       ("elsconvoldgsftkidx"    ).setBranchAlias("els_conv_old_gsftkidx"    );  //index of the GSF partner track, if thats where we find it
    //produces<vector<int>      >       ("elsconvolddelMissHits" ).setBranchAlias("els_conv_old_delMissHits" );  //Delta Missing Hits between the electron and partner track
    //produces<vector<int>      >       ("elsconvoldflag"        ).setBranchAlias("els_conv_old_flag"        );
    produces<vector<bool>     >       ("elsconvvtxflag"        ).setBranchAlias("els_conv_vtx_flag"        );
 
    // predefined 2012 ID decisions
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
    // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h?revision=1.5&view=markup
    //

    ///////////////////
    // Added for 53x //
    ///////////////////

    produces<vector<bool > >  ("elspassingMvaPreselection"  ).setBranchAlias("els_passingMvaPreselection"  );
    produces<vector<bool > >  ("elspassingPflowPreselection").setBranchAlias("els_passingPflowPreselection");
    produces<vector<float> >  ("elsr9"                      ).setBranchAlias("els_r9"                      );
    produces<vector<float> >  ("elssigmaIphiIphi"           ).setBranchAlias("els_sigmaIphiIphi"           );

    ///////////////////
    // Added for 7   //
    ///////////////////

    produces<vector<vector<int>   >   >       ("elspfcandidx"    ).setBranchAlias("els_PFCand_idx"    );
    produces<vector<float>            >       ("elsmass"         ).setBranchAlias("els_mass"          );

    //////////////////////
    // genMatch miniAOD //
    //////////////////////

  produces<vector<int>           >("elsmcpatMatchid"            	).setBranchAlias("els_mc_patMatch_id"          		); 
  produces<vector<LorentzVector> >("elsmcpatMatchp4"            	).setBranchAlias("els_mc_patMatch_p4"          		);
  produces<vector<float>         >("elsmcpatMatchdr"            	).setBranchAlias("els_mc_patMatch_dr"                  	);

  produces<vector<float>         >("elssigmaIPhiIPhifull5x5"   	).setBranchAlias("els_sigmaIPhiIPhi_full5x5"                  	);
  produces<vector<float>         >("elssigmaEtaEtafull5x5"     	).setBranchAlias("els_sigmaEtaEta_full5x5"                    	);
  produces<vector<float>         >("elssigmaIEtaIEtafull5x5"   	).setBranchAlias("els_sigmaIEtaIEta_full5x5"                  	);
  produces<vector<float>         >("elsr9full5x5"              	).setBranchAlias("els_r9_full5x5"                             	);
  produces<vector<float>         >("else1x5full5x5"            	).setBranchAlias("els_e1x5_full5x5"                           	);
  produces<vector<float>         >("else5x5full5x5"            	).setBranchAlias("els_e5x5_full5x5"                           	);
  produces<vector<float>         >("else2x5Maxfull5x5"         	).setBranchAlias("els_e2x5Max_full5x5"                        	);

  produces<vector<float>         >("elsminiIsouncor"       ).setBranchAlias("els_miniIso_uncor"                       	);
  produces<vector<float>         >("elsminiIsoch"       ).setBranchAlias("els_miniIso_ch"                       	);
  produces<vector<float>         >("elsminiIsonh"       ).setBranchAlias("els_miniIso_nh"                       	);
  produces<vector<float>         >("elsminiIsoem"       ).setBranchAlias("els_miniIso_em"                       	);
  produces<vector<float>         >("elsminiIsodb"       ).setBranchAlias("els_miniIso_db"                       	);

  produces<vector<float>         >("elsecalPFClusterIso"       ).setBranchAlias("els_ecalPFClusterIso"                       	);
  produces<vector<float>         >("elshcalPFClusterIso"       ).setBranchAlias("els_hcalPFClusterIso"                       	);


    // for matching to vertices using the "PFNoPileup" method
    // hint: it is just track vertex association 
    pfPileUpAlgo_ = new PFPileUpAlgo();
}

ElectronMaker::~ElectronMaker()
{
  if (pfPileUpAlgo_) delete pfPileUpAlgo_;
  if (clusterTools_) delete clusterTools_;
  if (mtsTransform_) delete mtsTransform_;
}

void  ElectronMaker::beginRun(const edm::Run&, const EventSetup& es) {
  
    ESHandle<TrackerGeometry>              trackerGeometryHandle;
    ESHandle<MagneticField>                magFieldHandle;
  
    es.get<TrackerDigiGeometryRecord>().get(trackerGeometryHandle);
    es.get<IdealMagneticFieldRecord>().get(magFieldHandle);
    mtsTransform_ = new MultiTrajectoryStateTransform(trackerGeometryHandle.product(), magFieldHandle.product());
}

void ElectronMaker::beginJob() {
}

void ElectronMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void ElectronMaker::produce(Event& iEvent, const EventSetup& iSetup) {

    // Define vectors to be filled
    auto_ptr<unsigned int>   evt_nels(new unsigned int) ;

    // ECAL related (superCluster) variables
    auto_ptr<vector<int> >   els_nSeed       (new vector<int>   );
    auto_ptr<vector<float> > els_etaSC       (new vector<float> );
    auto_ptr<vector<float> > els_phiSC       (new vector<float> );
    auto_ptr<vector<float> > els_eSC         (new vector<float> );
    auto_ptr<vector<float> > els_eSCRaw      (new vector<float> );
    auto_ptr<vector<float> > els_eSCPresh    (new vector<float> );
    auto_ptr<vector<int> >   els_fiduciality (new vector<int>   );
    auto_ptr<vector<int> >   els_type        (new vector<int>   );
//    auto_ptr<vector<int> >   els_scindex     (new vector<int>   ); 
    auto_ptr<vector<float> > els_etaSCwidth  (new vector<float> );
    auto_ptr<vector<float> > els_phiSCwidth  (new vector<float> );

    // uncertainties and corrections
    // somewhat complicated: see 
    // http://cms-service-sdtweb.web.cern.ch/cms-service-sdtweb/doxygen/CMSSW_3_1_2/doc/html/d5/d4b/GsfElectron_8h-source.html
    // note that if ecalEnergy == eSC depends on if further ecal corrections have been applied to the electron
    // after its construction
    auto_ptr<vector<float> > els_ecalEnergy            (new vector<float>);
    auto_ptr<vector<float> > els_ecalEnergyError       (new vector<float>);
    auto_ptr<vector<float> > els_trackMomentumError    (new vector<float>);
    auto_ptr<vector<float> > els_electronMomentumError (new vector<float>);
  
    // ID variables
    //
    //auto_ptr<vector<float> > els_lh                            (new vector<float> );
    auto_ptr<vector<float> > els_mva                           (new vector<float> );
    auto_ptr<vector<float> > els_dEtaIn                        (new vector<float> );
    auto_ptr<vector<float> > els_dEtaOut                       (new vector<float> );
    auto_ptr<vector<float> > els_dPhiIn                        (new vector<float> );
    auto_ptr<vector<float> > els_dPhiOut                       (new vector<float> );
    auto_ptr<vector<float> > els_dPhiInPhiOut                  (new vector<float> );
    auto_ptr<vector<float> > els_fbrem                         (new vector<float> );
    auto_ptr<vector<float> > els_eSeed                         (new vector<float> );
    auto_ptr<vector<float> > els_eOverPIn                      (new vector<float> );
    auto_ptr<vector<float> > els_eSeedOverPOut                 (new vector<float> );
    auto_ptr<vector<float> > els_eSeedOverPIn                  (new vector<float> );
    auto_ptr<vector<float> > els_eOverPOut                     (new vector<float> );
    auto_ptr<vector<float> > els_deltaEtaEleClusterTrackAtCalo (new vector<float> );
    auto_ptr<vector<float> > els_deltaPhiEleClusterTrackAtCalo (new vector<float> );
    auto_ptr<vector<bool > > els_isGsfCtfScPixChargeConsistent (new vector<bool > );
                             
    auto_ptr<vector<float> > els_hOverE                        (new vector<float> );
    auto_ptr<vector<float> > els_hOverEBC                      (new vector<float> );
    auto_ptr<vector<float> > els_hcalDepth1OverEcal            (new vector<float> );
    auto_ptr<vector<float> > els_hcalDepth2OverEcal            (new vector<float> );
                             
    //auto_ptr<vector<float> > els_sigmaPhiPhi                   (new vector<float> );
    auto_ptr<vector<float> > els_sigmaIPhiIPhi                 (new vector<float> );
    //auto_ptr<vector<float> > els_sigmaIEtaIPhi                 (new vector<float> );
    auto_ptr<vector<float> > els_sigmaEtaEta                   (new vector<float> );
    auto_ptr<vector<float> > els_sigmaIEtaIEta                 (new vector<float> );
    //auto_ptr<vector<float> > els_sigmaIPhiIPhiSC               (new vector<float> );
    //auto_ptr<vector<float> > els_sigmaIEtaIEtaSC               (new vector<float> );
                             
    auto_ptr<vector<float> > els_e2x5Max                       (new vector<float> );
    auto_ptr<vector<float> > els_e1x5                          (new vector<float> );
    auto_ptr<vector<float> > els_e5x5                          (new vector<float> );
//    auto_ptr<vector<float> > els_e3x3                          (new vector<float> );
//    auto_ptr<vector<float> > els_eMax                          (new vector<float> );

    // predefined ID decisions
    //
    auto_ptr<vector<int> > els_class    (new vector<int>);
    auto_ptr<vector<int> > els_category (new vector<int>);

    auto_ptr<vector<int> > passVetoId     (new vector<int>);
    auto_ptr<vector<int> > passLooseId    (new vector<int>);
    auto_ptr<vector<int> > passMediumId   (new vector<int>);
    auto_ptr<vector<int> > passTightId    (new vector<int>);

    // isolation variables
    //
    auto_ptr<vector<float> > els_tkIso                  (new vector<float> );
    auto_ptr<vector<float> > els_ecalIso                (new vector<float> );
    auto_ptr<vector<float> > els_hcalIso                (new vector<float> );
    auto_ptr<vector<float> > els_hcalDepth1TowerSumEt   (new vector<float> );
    auto_ptr<vector<float> > els_hcalDepth2TowerSumEt   (new vector<float> );
                             
    auto_ptr<vector<float> > els_tkIso04                (new vector<float> );
    auto_ptr<vector<float> > els_ecalIso04              (new vector<float> );
    auto_ptr<vector<float> > els_hcalIso04              (new vector<float> );
    auto_ptr<vector<float> > els_hcalDepth1TowerSumEt04 (new vector<float> );
    auto_ptr<vector<float> > els_hcalDepth2TowerSumEt04 (new vector<float> );

//    auto_ptr<vector<float> > els_iso03_pf               (new vector<float> );
//    auto_ptr<vector<float> > els_iso04_pf               (new vector<float> );

//    auto_ptr<vector<float> > els_iso03_pf_ch            (new vector<float> );
//    auto_ptr<vector<float> > els_iso03_pf_gamma05       (new vector<float> );
//    auto_ptr<vector<float> > els_iso03_pf_nhad05        (new vector<float> );
//    auto_ptr<vector<float> > els_iso04_pf_ch            (new vector<float> );
//    auto_ptr<vector<float> > els_iso04_pf_gamma05       (new vector<float> );
//    auto_ptr<vector<float> > els_iso04_pf_nhad05        (new vector<float> );

    auto_ptr<vector<float> > els_pfChargedHadronIso     (new vector<float> );
    auto_ptr<vector<float> > els_pfNeutralHadronIso     (new vector<float> );
    auto_ptr<vector<float> > els_pfPhotonIso            (new vector<float> );
    auto_ptr<vector<float> > els_pfPUIso                (new vector<float> );

    // track variables
    //
    auto_ptr<vector<int> >   els_charge     (new vector<int>   );
    auto_ptr<vector<int> >   els_trk_charge (new vector<int>   );
    auto_ptr<vector<int> >   els_sccharge   (new vector<int>   );
    auto_ptr<vector<float> > els_chi2       (new vector<float> );
    auto_ptr<vector<float> > els_ndof       (new vector<float> );
    auto_ptr<vector<float> > els_ckf_chi2   (new vector<float> );
    auto_ptr<vector<float> > els_ckf_ndof   (new vector<float> );
    auto_ptr<vector<int> >   els_validHits  (new vector<int>   );
    auto_ptr<vector<int> >   els_lostHits   (new vector<int>   );
    auto_ptr<vector<float> > els_d0         (new vector<float> );
    auto_ptr<vector<float> > els_z0         (new vector<float> );
    auto_ptr<vector<float> > els_dxyPV      (new vector<float> );
    auto_ptr<vector<float> > els_dzPV       (new vector<float> );
    auto_ptr<vector<float> > els_d0Err      (new vector<float> );
    auto_ptr<vector<float> > els_z0Err      (new vector<float> );
    auto_ptr<vector<float> > els_d0corr     (new vector<float> );
    auto_ptr<vector<float> > els_d0corrPhi  (new vector<float> );
    auto_ptr<vector<float> > els_d0phiCov   (new vector<float> );
    auto_ptr<vector<float> > els_z0corr     (new vector<float> );
    auto_ptr<vector<float> > els_ptErr      (new vector<float> );
    auto_ptr<vector<float> > els_ptErrGsf   (new vector<float> );
    auto_ptr<vector<float> > els_etaErr     (new vector<float> );
    auto_ptr<vector<float> > els_phiErr     (new vector<float> );
    auto_ptr<vector<float> > els_ip3d       (new vector<float> );
    auto_ptr<vector<float> > els_ip3derr    (new vector<float> );
    auto_ptr<vector<float> > els_ip2d       (new vector<float> );
    auto_ptr<vector<float> > els_ip2derr    (new vector<float> );
    auto_ptr<vector<float> > els_bs3d       (new vector<float> );
    auto_ptr<vector<float> > els_bs3derr    (new vector<float> );
    auto_ptr<vector<float> > els_bs2d       (new vector<float> );
    auto_ptr<vector<float> > els_bs2derr    (new vector<float> );
  
    // LorentzVectors
    //
    auto_ptr<vector<LorentzVector> > els_p4     (new vector<LorentzVector>);
    auto_ptr<vector<LorentzVector> > els_trk_p4 (new vector<LorentzVector>);
    auto_ptr<vector<LorentzVector> > els_p4In   (new vector<LorentzVector>);
    auto_ptr<vector<LorentzVector> > els_p4Out  (new vector<LorentzVector>);

    // Vertex
    //
    auto_ptr<vector<LorentzVector> > els_vertex_p4 (new vector<LorentzVector>);
    auto_ptr<vector<LorentzVector> > els_trk_vertex_p4 (new vector<LorentzVector>);

    //HitPattern information
    //
    auto_ptr<vector<LorentzVector> >          els_inner_position       (new vector<LorentzVector> );
    auto_ptr<vector<LorentzVector> >          els_outer_position       (new vector<LorentzVector> );
    auto_ptr<vector<int> >                    els_valid_pixelhits      (new vector<int>           ); 
    auto_ptr<vector<int> >                    els_lost_pixelhits       (new vector<int>           ); 
    auto_ptr<vector<int> >                    els_nlayers              (new vector<int>           ); 
    auto_ptr<vector<int> >                    els_nlayers3D            (new vector<int>           ); 
    auto_ptr<vector<int> >                    els_nlayersLost          (new vector<int>           ); 
    //auto_ptr<vector<int> >                    els_layer1_sizerphi      (new vector<int>           ); 
    //auto_ptr<vector<int> >                    els_layer1_sizerz        (new vector<int>           ); 
    //auto_ptr<vector<float> >                  els_layer1_charge        (new vector<float>         );
    //auto_ptr<vector<int> >                    els_layer1_det           (new vector<int>           );
    //auto_ptr<vector<int> >                    els_layer1_layer         (new vector<int>           );
    auto_ptr<vector<int> >                    els_exp_innerlayers      (new vector<int>           ); 
    auto_ptr<vector<int> >                    els_exp_outerlayers      (new vector<int>           ); 
    auto_ptr<vector<int> >                    els_ckf_laywithmeas      (new vector<int>           );
    auto_ptr<vector<int> >                    els_ckf_charge           (new vector<int>           );

    auto_ptr<vector<float> >                  els_trkshFrac            (new vector<float>         );
    auto_ptr<vector<float> >                  els_trkdr                (new vector<float>         );

    //conversions
    //auto_ptr<vector<vector<float> > >         els_convs_dist           (new vector<vector<float> > );
    //auto_ptr<vector<vector<float> > >         els_convs_dcot           (new vector<vector<float> > );
    //auto_ptr<vector<vector<float> > >         els_convs_radius         (new vector<vector<float> > );
    //auto_ptr<vector<vector<int> > >           els_convs_tkidx          (new vector<vector<int> >   );
    //auto_ptr<vector<vector<int> > >           els_convs_gsftkidx       (new vector<vector<int> >   );
    //auto_ptr<vector<vector<int> > >           els_convs_delMissHits    (new vector<vector<int> >   );
    //auto_ptr<vector<vector<int> > >           els_convs_flag           (new vector<vector<int> >   );  
    //  
    //auto_ptr<vector<LorentzVector> >          els_conv_pos_p4          (new vector<LorentzVector>  );
    //auto_ptr<vector<float> >                  els_conv_dist            (new vector<float>          );
    //auto_ptr<vector<float> >                  els_conv_dcot            (new vector<float>          );
    //auto_ptr<vector<float> >                  els_conv_radius          (new vector<float>          );
    //auto_ptr<vector<int>   >                  els_conv_tkidx           (new vector<int>            );
    //auto_ptr<vector<int>   >                  els_conv_gsftkidx        (new vector<int>            );
    //auto_ptr<vector<int>   >                  els_conv_delMissHits     (new vector<int>            );
    //auto_ptr<vector<int>   >                  els_conv_flag            (new vector<int>            );

    //auto_ptr<vector<float> >                  els_conv_old_dist        (new vector<float>          );
    //auto_ptr<vector<float> >                  els_conv_old_dcot        (new vector<float>          );
    //auto_ptr<vector<float> >                  els_conv_old_radius      (new vector<float>          );
    //auto_ptr<vector<int> >                    els_conv_old_tkidx       (new vector<int>            );
    //auto_ptr<vector<int> >                    els_conv_old_gsftkidx    (new vector<int>            );
    //auto_ptr<vector<int> >                    els_conv_old_delMissHits (new vector<int>            );
    //auto_ptr<vector<int> >                    els_conv_old_flag        (new vector<int>            );

    auto_ptr<vector<bool> >                   els_conv_vtx_flag        (new vector<bool>           );

    ///////////////////
    // Added for 53x //
    ///////////////////

    auto_ptr<vector<bool > >  els_passingMvaPreselection   ( new vector<bool>  );
    auto_ptr<vector<bool > >  els_passingPflowPreselection ( new vector<bool>  );
    auto_ptr<vector<float> >  els_r9                       ( new vector<float> );
    auto_ptr<vector<float> >  els_sigmaIphiIphi            ( new vector<float> );

    ///////////////////
    // Added for 7   //
    ///////////////////

    auto_ptr<vector<vector<int> > >           els_PFCand_idx       (new vector<vector<int> >   );
    auto_ptr<vector<float> >                  els_mass             (new vector<float>          );

    //////////////////////
    // Added miniAOD    //
    //////////////////////
  auto_ptr<vector<int>           >       els_mc_patMatch_id          (new vector<int>          );
  auto_ptr<vector<LorentzVector> >       els_mc_patMatch_p4          (new vector<LorentzVector>);
  auto_ptr<vector<float>         >       els_mc_patMatch_dr          (new vector<float>        );

  auto_ptr<vector<float>   >       els_sigmaIPhiIPhi_full5x5             (new vector<float>        );
  auto_ptr<vector<float>   >       els_sigmaEtaEta_full5x5               (new vector<float>        );
  auto_ptr<vector<float>   >       els_sigmaIEtaIEta_full5x5             (new vector<float>        );
  auto_ptr<vector<float>   >       els_r9_full5x5                        (new vector<float>        );
  auto_ptr<vector<float>   >       els_e1x5_full5x5                      (new vector<float>        );
  auto_ptr<vector<float>   >       els_e5x5_full5x5                      (new vector<float>        );
  auto_ptr<vector<float>   >       els_e2x5Max_full5x5                   (new vector<float>        );

  auto_ptr<vector<float>   >       els_miniIso_uncor                  (new vector<float>        );  	
  auto_ptr<vector<float>   >       els_miniIso_ch                  (new vector<float>        );  	
  auto_ptr<vector<float>   >       els_miniIso_nh                  (new vector<float>        );  	
  auto_ptr<vector<float>   >       els_miniIso_em                  (new vector<float>        );  	
  auto_ptr<vector<float>   >       els_miniIso_db                  (new vector<float>        );  	

  auto_ptr<vector<float>   >       els_ecalPFClusterIso                (new vector<float>        );  	
  auto_ptr<vector<float>   >       els_hcalPFClusterIso                (new vector<float>        );  	


    // --- Get Input Collections --- //

    ////////////////
    // Get Tracks //
    ////////////////
   
//    Handle<TrackCollection> tracks_h;
//    iEvent.getByLabel(trksInputTag_, tracks_h);

  
    ////////////////
    // GSF Tracks //
    ////////////////

//    Handle<GsfTrackCollection> gsftracks_h;
//    iEvent.getByLabel(gsftracksInputTag_, gsftracks_h);


    /////////////
    // B Field //
    /////////////

    Handle<float> evt_bField_h;
    iEvent.getByLabel("eventMaker", "evtbField", evt_bField_h);
    float evt_bField = *evt_bField_h.product();
    if ( evt_bField == 1234567 ) ; // Avoid "unused variable" error while the function using this variable is inactive
    
    ///////////////
    // Electrons //
    ///////////////

    Handle<View<pat::Electron> > els_h;
    iEvent.getByLabel(electronsInputTag_, els_h);
    View<pat::Electron> gsfElColl = *(els_h.product());



//    Handle<GsfElectronCollection> els_coll_h;
//    iEvent.getByLabel(electronsInputTag_, els_coll_h);    

    //////////////
    // PF Cands //
    //////////////

     iEvent.getByLabel(pfCandsInputTag, packPfCand_h);
     pfCandidates  = packPfCand_h.product();

    /////////////////////////
    // External Isolations //
    /////////////////////////
    // edm::Handle< edm::ValueMap<double> > pfIsoCharged03_h;
    // iEvent.getByLabel(pfIsoCharged03InputTag, pfIsoCharged03_h);
    // edm::Handle< edm::ValueMap<double> > pfIsoGamma03_h;
    // iEvent.getByLabel(pfIsoGamma03InputTag, pfIsoGamma03_h);
    // edm::Handle< edm::ValueMap<double> > pfIsoNeutral03_h;
    // iEvent.getByLabel(pfIsoNeutral03InputTag, pfIsoNeutral03_h);

    // edm::Handle< edm::ValueMap<double> > pfIsoCharged04_h;
    // iEvent.getByLabel(pfIsoCharged04InputTag, pfIsoCharged04_h);
    // edm::Handle< edm::ValueMap<double> > pfIsoGamma04_h;
    // iEvent.getByLabel(pfIsoGamma04InputTag, pfIsoGamma04_h);
    // edm::Handle< edm::ValueMap<double> > pfIsoNeutral04_h;
    // iEvent.getByLabel(pfIsoNeutral04InputTag, pfIsoNeutral04_h);

  
    ////////////
    // Vertex //
    ////////////

    iEvent.getByLabel(vtxInputTag, vertexHandle);


    ///////////////////////////
    // TransientTrackBuilder //
    ///////////////////////////
//    ESHandle<TransientTrackBuilder> theTTBuilder;
//    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);


    ////////////////////////////////////////////////
    // Get tools to get cluster shape information //
    ////////////////////////////////////////////////

//    if ( clusterTools_ ) delete clusterTools_;
//    clusterTools_ = new EcalClusterLazyTools( iEvent, iSetup, InputTag("reducedEcalRecHitsEB"), InputTag("reducedEcalRecHitsEE") );


    //////////////
    // Beamspot //
    //////////////

    InputTag beamSpot_tag(beamSpotInputTag_.label(),"evtbsp4");
    Handle<LorentzVector> beamSpotH;
    iEvent.getByLabel(beamSpot_tag, beamSpotH);
    const Point beamSpot = beamSpotH.isValid() ? Point(beamSpotH->x(), beamSpotH->y(), beamSpotH->z()) : Point(0,0,0);

    Handle<reco::BeamSpot> beamspot_h;
    iEvent.getByLabel(beamSpot_tag_, beamspot_h);
    const reco::BeamSpot &beamSpotreco = *(beamspot_h.product()); 
    if ( beamSpotreco.x0() == 1234567 ) ; // Avoid "unused variable" error while the function using this variable is inactive


    ///////////////////////
    // rho for isolation //
    ///////////////////////

//    edm::Handle<float> rhoIso_h;
//    iEvent.getByLabel(rhoInputTag_, rhoIso_h);
//    float rhoIso = *(rhoIso_h.product());

    /////////////////////////////////////////////////////////////
    // Get the value maps for the Egamma electron ID decisions //
    /////////////////////////////////////////////////////////////

//    const ValueMap<float>&  eidLHMap = getValueMap<float>(iEvent, eidLHTag_);

  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);


    //////////////////////////
    // get cms2scsseeddetid //
    //////////////////////////

//    InputTag cms2scsseeddetid_tag(cms2scsseeddetidInputTag_.label(),"scsdetIdSeed");
//    Handle<vector<int> > cms2scsseeddetid_h;
//    iEvent.getByLabel(cms2scsseeddetid_tag, cms2scsseeddetid_h); 
//    const vector<int> *cms2scsseeddetid = cms2scsseeddetid_h.product();

    //////////////////////////////
    // Get the ele<->PFCand map //
    //////////////////////////////

//    edm::Handle<edm::ValueMap<std::vector<reco::PFCandidateRef > > > eleToParticleBasedIsoMapHandle;
//    InputTag particleBase(string("particleBasedIsolation"),string("gedGsfElectrons"));  
//    iEvent.getByLabel(particleBase, eleToParticleBasedIsoMapHandle);    
//    edm::ValueMap<std::vector<reco::PFCandidateRef > >   eleToParticleBasedIsoMap =  *(eleToParticleBasedIsoMapHandle.product());
    
    // --- Fill --- //

    /////////////////////////
    // Loop Over Electrons //
    /////////////////////////

    *evt_nels       = els_h->size();
    double mass     = 0.000510998918;
    size_t elsIndex = 0;
    for( View<pat::Electron>::const_iterator el = els_h->begin(); el != els_h->end(); el++, elsIndex++ ) {


        ////////////////
        // References //
        ////////////////
      const GsfTrackRef            el_track         = el->gsfTrack(); // Embedded GSF Track for miniAOD
      const RefToBase<pat::Electron> gsfElRef         = els_h->refAt(elsIndex);    
      const TrackRef               ctfTkRef         = el->closestCtfTrackRef(); // Embedded CTF Track for miniAOD 

/*
        const Track*                 el_track         = (const Track*)(el->gsfTrack().get());
        const RefToBase<pat::Electron> gsfElRef         = els_h->refAt(elsIndex);    

        //const TrackRef               ctfTkRef         = el->closestCtfTrackRef();
        const GsfTrackRef            gsfTkRef         = el->gsfTrack();

*/

        ////////////
        // Vertex //
        ////////////
        const VertexCollection*      vertexCollection = vertexHandle.product();
        VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
        int firstGoodVertexIdx = 0;
        for (VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++firstGoodVertexIdx) {
	  // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
	  // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
	  if (  /*!vtx->isFake() &&*/ !(vtx->chi2()==0 && vtx->ndof()==0) &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
                firstGoodVertex = vtx;
                break;
	  }
        }

        //////////////////////
        // Fiduciality Mask //
        //////////////////////

        int fiducialityMask = 0;  // the enum is in interface/EgammaFiduciality.h
        if ( el->isEB()        ) fiducialityMask |= 1 << ISEB;
        if ( el->isEBEEGap()   ) fiducialityMask |= 1 << ISEBEEGAP;
        if ( el->isEE()        ) fiducialityMask |= 1 << ISEE;
        if ( el->isEEGap()     ) fiducialityMask |= 1 << ISEEGAP;
        if ( el->isEBEtaGap()  ) fiducialityMask |= 1 << ISEBETAGAP;
        if ( el->isEBPhiGap()  ) fiducialityMask |= 1 << ISEBPHIGAP;
        if ( el->isEEDeeGap()  ) fiducialityMask |= 1 << ISEEDEEGAP;
        if ( el->isEERingGap() ) fiducialityMask |= 1 << ISEERINGGAP;
        if ( el->isGap()       ) fiducialityMask |= 1 << ISGAP;

  
        ///////////////////////////
        // Corrections & Seeding //
        ///////////////////////////

        int electronTypeMask = 0;
        if ( el->isEcalEnergyCorrected()        ) electronTypeMask |= 1 << ISECALENERGYCORRECTED;
        if ( el->trackerDrivenSeed()            ) electronTypeMask |= 1 << ISTRACKERDRIVEN;
        if ( el->ecalDrivenSeed()               ) electronTypeMask |= 1 << ISECALDRIVEN;
        if ( el->passingCutBasedPreselection()  ) electronTypeMask |= 1 << ISCUTPRESELECTED;
        if ( el->passingMvaPreselection()       ) electronTypeMask |= 1 << ISMVAPRESELECTED;
        //if ( el->isMomentumCorrected() ) electronTypeMask |= 1 << ISMOMENTUMCORRECTED; // Depricated in CMSSW_4_2x ( DataFormats/EgammaCandidates/interface/GsfElectron.h )

        /////////////////////
        // Lorentz Vectors //
        /////////////////////

        LorentzVector    p4In;
        LorentzVector    p4Out;
        LorentzVector    trk_p4( el_track->px(), el_track->py(), el_track->pz(), el_track->p() );
        math::XYZVectorF p3In  = el->trackMomentumAtVtx();
        math::XYZVectorF p3Out = el->trackMomentumOut();
        p4In.SetXYZT (   p3In.x() , p3In.y() , p3In.z() , sqrt( mass*mass + p3In.R() *p3In.R()  ) );
        p4Out.SetXYZT(   p3Out.x(), p3Out.y(), p3Out.z(), sqrt( mass*mass + p3Out.R()*p3Out.R() ) );


        ///////////////////
        // Predifined ID //
        ///////////////////

        els_class              ->push_back( el->classification()  ); // this is the old pTDR classification
        els_category           ->push_back( classify(gsfElRef)    ); // this is the sani classification

        const Ptr<pat::Electron> elPtr(els_h, el - els_h->begin() );
  
        passVetoId  ->push_back( (*veto_id_decisions)[ elPtr ] );
        passLooseId ->push_back( (*loose_id_decisions)[ elPtr ] );
        passMediumId->push_back( (*medium_id_decisions)[ elPtr ] );
        passTightId ->push_back( (*tight_id_decisions)[ elPtr ] );


        //////////////
        // Electron //
        //////////////

        els_fiduciality        ->push_back( fiducialityMask                                 );
        els_type               ->push_back( electronTypeMask                                );
        //els_ecalEnergy         ->push_back( el->ecalEnergy()                                );  // energy corrections and uncertainties
        //els_ecalEnergyError    ->push_back( el->ecalEnergyError()                           );
        els_ecalEnergy         ->push_back( el->correctedEcalEnergy()                       );  // energy corrections and uncertainties
        els_ecalEnergyError    ->push_back( el->correctedEcalEnergyError()                  );
        els_trackMomentumError ->push_back( el->trackMomentumError()                        );
        els_p4                 ->push_back( LorentzVector( el->p4() )                       );
        els_trk_p4             ->push_back( trk_p4                                          );
        els_p4In               ->push_back( p4In                                            );
        els_p4Out              ->push_back( p4Out                                           );
        els_vertex_p4          ->push_back( LorentzVector(el->vx(), el->vy(), el->vz(), 0.) );
        els_trk_vertex_p4      ->push_back( LorentzVector(el_track->vx(), el_track->vy(), el_track->vz(), 0.) );


        ///////////////
        // Isolation //
        ///////////////

        els_ecalIso               ->push_back( el->dr03EcalRecHitSumEt()                  );
        els_hcalIso               ->push_back( el->dr03HcalTowerSumEt()                   );
        els_hcalDepth1TowerSumEt  ->push_back( el->dr03HcalDepth1TowerSumEt()             );
        els_hcalDepth2TowerSumEt  ->push_back( el->dr03HcalDepth2TowerSumEt()             );
        els_tkIso                 ->push_back( el->dr03TkSumPt()                          );

        els_ecalIso04             ->push_back( el->dr04EcalRecHitSumEt()                  );
        els_hcalIso04             ->push_back( el->dr04HcalTowerSumEt()                   );
        els_hcalDepth1TowerSumEt04->push_back( el->dr04HcalDepth1TowerSumEt()             );
        els_hcalDepth2TowerSumEt04->push_back( el->dr04HcalDepth2TowerSumEt()             );
        els_tkIso04               ->push_back( el->dr04TkSumPt()                          );

    
        //////////////////
        // PF Isolation //
        //////////////////

        GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();

        els_pfChargedHadronIso -> push_back( pfIso.sumChargedHadronPt );
        els_pfNeutralHadronIso -> push_back( pfIso.sumNeutralHadronEt );
        els_pfPhotonIso        -> push_back( pfIso.sumPhotonEt        );
        els_pfPUIso            -> push_back( pfIso.sumPUPt            );


        //////////////////
        // Supercluster //
        //////////////////

        els_etaSC         ->push_back( el->superCluster()->eta()             );
        els_phiSC         ->push_back( el->superCluster()->phi()             );
	els_eSC           ->push_back( el->superCluster()->energy()          );
        els_eSCRaw        ->push_back( el->superCluster()->rawEnergy()       );
        els_eSCPresh      ->push_back( el->superCluster()->preshowerEnergy() );
        els_nSeed         ->push_back( el->basicClustersSize() - 1           );
        els_e1x5          ->push_back( el->e1x5()                            );
        els_e5x5          ->push_back( el->e5x5()                            );
        els_e2x5Max       ->push_back( el->e2x5Max()                         );
        els_sigmaEtaEta   ->push_back( el->sigmaEtaEta()                     );
        els_sigmaIEtaIEta ->push_back( el->sigmaIetaIeta()                   );
        els_etaSCwidth    ->push_back( el->superCluster()->etaWidth()        );
        els_phiSCwidth    ->push_back( el->superCluster()->phiWidth()        );
	
	// We used to make these using the cluster tools, but now we can take them directly from RECO electron
	els_sigmaIPhiIPhi         ->push_back( el->sigmaIphiIphi()           );

	// Take these directly from the PAT electron of the miniAOD
	els_sigmaIPhiIPhi_full5x5  ->push_back( el->full5x5_sigmaIphiIphi()   ); 
	els_sigmaEtaEta_full5x5    ->push_back( el->full5x5_sigmaEtaEta()     );
	els_sigmaIEtaIEta_full5x5  ->push_back( el->full5x5_sigmaIetaIeta()   );
	els_r9_full5x5             ->push_back( el->full5x5_r9()              );
	els_e1x5_full5x5           ->push_back( el->full5x5_e1x5()            );
	els_e5x5_full5x5           ->push_back( el->full5x5_e5x5()            );  
	els_e2x5Max_full5x5        ->push_back( el->full5x5_e2x5Max()         );



        ///////////////////////////////////////////////////////
        // Get cluster info that is not stored in the object //
        ///////////////////////////////////////////////////////

        if( el->superCluster()->seed().isAvailable() ) { 
//
//            //
//            const BasicCluster&  clRef              = *(el->superCluster()->seed());
//            const vector<float>& covs               = clusterTools_->covariances(clRef);                         // get the covariances computed in 5x5 around the seed
//            const vector<float>& lcovs              = clusterTools_->localCovariances(clRef);                    // get the local covariances computed in a 5x5 around the seed
//            const vector<float>  localCovariancesSC = clusterTools_->scLocalCovariances(*(el->superCluster()));  // get the local covariances computed using all crystals in the SC
//
//            //
            els_eSeed           ->push_back( el->superCluster()->seed()->energy()     );
//            els_sigmaPhiPhi     ->push_back( isfinite(covs[2])               ? covs[2] > 0                ? sqrt(covs[2])  : -1 * sqrt(-1 * covs[2])                              : -9999. );
//get from RECO            els_sigmaIPhiIPhi   ->push_back( isfinite(lcovs[2])              ? lcovs[2] > 0               ? sqrt(lcovs[2]) : -1 * sqrt(-1 * lcovs[2])                             : -9999. );
//            els_sigmaIEtaIPhi   ->push_back( isfinite(lcovs[1])              ? lcovs[1] > 0               ? sqrt(lcovs[1]) : -1 * sqrt(-1 * lcovs[1])                             : -9999. );
//get from RECO            els_sigmaIEtaIEtaSC ->push_back( isfinite(localCovariancesSC[0]) ? localCovariancesSC[0] > 0  ? sqrt(localCovariancesSC[0])   : -1 * sqrt(-1 * localCovariancesSC[0]) : -9999. );
//            els_sigmaIPhiIPhiSC ->push_back( isfinite(localCovariancesSC[2]) ? localCovariancesSC[2] > 0  ? sqrt(localCovariancesSC[2])   : -1 * sqrt(-1 * localCovariancesSC[2]) : -9999. );
//
//            //
//            els_e3x3            ->push_back( clusterTools_->e3x3(clRef) );
//            els_eMax            ->push_back( clusterTools_->eMax(clRef) );
        } 
        else {

            //
            els_eSeed           ->push_back(-9999.);
//            els_sigmaPhiPhi     ->push_back(-9999.);
            els_sigmaIPhiIPhi   ->push_back(-9999.);
//            els_sigmaIEtaIPhi   ->push_back(-9999.);
//            els_sigmaIEtaIEtaSC ->push_back(-9999.);
//            els_sigmaIPhiIPhiSC ->push_back(-9999.);
//
//            //
//            els_e3x3            ->push_back(-9999.);
//            els_eMax            ->push_back(-9999.);
//
        } //
 
   
//        /////////////////////////
//        // Super Cluster Index //
//        /////////////////////////
//
//        // get cms2scsseeddetid--sorry for junk from photons...
//        bool foundseed = false;
//        for( unsigned int i=0; i<cms2scsseeddetid->size(); i++ ) {      
//
//            if( cms2scsseeddetid->at(i) == -9999            ) continue;      
//            if( !(el->superCluster()->seed().isAvailable()) ) continue;
//
//            if( uint32_t( cms2scsseeddetid->at(i) ) == el->superCluster()->seed()->seed() ) {
//                foundseed = true;
//                els_scindex->push_back( i );
//                break;
//            }
//
//        }
//        if( !foundseed ) {
//            els_scindex->push_back( -1 );
//            // this is understood: the photon can have energy significantly higher than SC for whatever reason.
//            // cout << "No seed found. seed id: " << int(photon->superCluster()->seed()->seed()) << "  photon et: " << photon->et() << "  sc et: " << photon->superCluster()->energy()/cosh(photon->superCluster()->eta()) << endl;
//        }


        ////////
        // ID //
        ////////

        double phi_pin  = ( el->caloPosition().phi() - el->deltaPhiSuperClusterTrackAtVtx() );
        double phi_pout = ( el->superCluster()->seed().isAvailable() ? ( el->superCluster()->seed()->position().phi() - el->deltaPhiSeedClusterTrackAtCalo() ) :  -9999. );

        //els_hOverE                        ->push_back( el->hadronicOverEm()                 );
        els_hOverE                        ->push_back( el->hcalOverEcal()                   );
        els_hOverEBC                      ->push_back( el->hcalOverEcalBc()                 );
        els_hcalDepth1OverEcal            ->push_back( el->hcalDepth1OverEcal()             );
        els_hcalDepth2OverEcal            ->push_back( el->hcalDepth2OverEcal()             );
        els_eOverPIn                      ->push_back( el->eSuperClusterOverP()             );
        els_eSeedOverPOut                 ->push_back( el->eSeedClusterOverPout()           );
        els_eSeedOverPIn                  ->push_back( el->eSeedClusterOverP()              );
        els_eOverPOut                     ->push_back( el->eEleClusterOverPout()            );
        els_fbrem                         ->push_back( el->fbrem()                          );
	//        els_lh                            ->push_back( eidLHMap[gsfElRef]                   );
        //els_mva                           ->push_back( el->mva()                            );
        els_mva                           ->push_back( el->mvaOutput().mva_Isolated                  );

        els_dEtaIn                        ->push_back( el->deltaEtaSuperClusterTrackAtVtx() );
        els_dEtaOut                       ->push_back( el->deltaEtaSeedClusterTrackAtCalo() );
        els_deltaEtaEleClusterTrackAtCalo ->push_back( el->deltaEtaEleClusterTrackAtCalo()  );
        els_dPhiIn                        ->push_back( el->deltaPhiSuperClusterTrackAtVtx() );
        els_dPhiInPhiOut                  ->push_back( phi_pin - phi_pout                   );
        els_dPhiOut                       ->push_back( el->deltaPhiSeedClusterTrackAtCalo() );
        els_deltaPhiEleClusterTrackAtCalo ->push_back( el->deltaPhiEleClusterTrackAtCalo()  );
	els_isGsfCtfScPixChargeConsistent ->push_back( el->isGsfCtfScPixChargeConsistent()  );



        
        ////////////
        // Tracks //
        ////////////

        float pt       = el_track->pt();
        float p        = el_track->p();
        float q        = el_track->charge();
        float pz       = el_track->pz();
        float trkpterr = (el_track->charge()!=0) ? sqrt(pt*pt*p*p/pow(q, 2)*(el_track->covariance(0,0))+2*pt*p/q*pz*(el_track->covariance(0,1))+ pz*pz*(el_track->covariance(1,1) ) ) : -9999.;
            
        els_chi2                  ->push_back( el_track->chi2()                          );
        els_ndof                  ->push_back( el_track->ndof()                          );
        els_d0Err                 ->push_back( el_track->d0Error()                       );
        els_z0Err                 ->push_back( el_track->dzError()                       );
        els_ptErr                 ->push_back( trkpterr                                  );
        els_ptErrGsf              ->push_back( el_track->ptError()                       );
        els_etaErr                ->push_back( el_track->etaError()                      );
        els_phiErr                ->push_back( el_track->phiError()                      );  
        els_validHits             ->push_back( el_track->numberOfValidHits()             );
        els_lostHits              ->push_back( el_track->numberOfLostHits()              );
        els_charge                ->push_back( el->charge()                              );
        els_trk_charge            ->push_back( el_track->charge()                        );
        els_sccharge              ->push_back( el->scPixCharge()                         );
        els_d0                    ->push_back( el_track->d0()                            );
        els_z0                    ->push_back( el_track->dz()                            );
        els_d0corr                ->push_back( -1*(el_track->dxy(beamSpot))              );
        els_z0corr                ->push_back( el_track->dz(beamSpot)                    );
	float d0corr = -1*(el_track->dxy(beamSpot));
	float corrd0phi = atan2( (-1 * d0corr * sin( el_track->phi() )), d0corr * cos( el_track->phi() ) );  
        els_d0corrPhi             ->push_back( corrd0phi                                 );
        els_d0phiCov              ->push_back( -1.* el_track->covariance(TrackBase::i_phi, TrackBase::i_dxy));
	if (firstGoodVertex!=vertexCollection->end()) {
	  els_dxyPV                 ->push_back( el_track->dxy( firstGoodVertex->position() )  );
	  els_dzPV                  ->push_back( el_track->dz(  firstGoodVertex->position() )  );
	}
	else {
	  els_dxyPV ->push_back( -999. );
	  els_dzPV  ->push_back( -999. );
	}

        /////////
        // CTF //
        /////////

        if( ctfTkRef.isNonnull() ) {
            //els_trkshFrac -> push_back( static_cast<float>( el->shFracInnerHits() )                                  );
            els_trkshFrac -> push_back( static_cast<float>( el->ctfGsfOverlap() )                                    );
            els_trkdr     -> push_back( deltaR( el_track->eta(), el_track->phi(), ctfTkRef->eta(), ctfTkRef->phi() ) );
            els_ckf_chi2  -> push_back( ctfTkRef->chi2() );
            els_ckf_ndof  -> push_back( ctfTkRef->ndof() );
            els_ckf_laywithmeas -> push_back( ctfTkRef->hitPattern().trackerLayersWithMeasurement() );
	    els_ckf_charge-> push_back( ctfTkRef->charge() );
        } 
        else {
            els_trkshFrac -> push_back(-9999.);
            els_trkdr     -> push_back(-9999.);
            els_ckf_chi2  -> push_back(-9999.);
            els_ckf_ndof  -> push_back(-9999.);
            els_ckf_laywithmeas -> push_back(-9999.);
	    els_ckf_charge-> push_back( 9999 );
        }

        
//        ////////////////////
//        // Regular Vertex //
//        ////////////////////        
//        TransientTrack tt = theTTBuilder->build(el->gsfTrack());
//    
//        if ( firstGoodVertex!=vertexCollection->end() ) {
//            Measurement1D ip3D_regular = IPTools::absoluteImpactParameter3D(tt, *firstGoodVertex).second;
//            //
//            els_ip3d      -> push_back( ip3D_regular.value() );
//            els_ip3derr   -> push_back( ip3D_regular.error() );
//        } else {
//            //
//            els_ip3d      -> push_back( -999. );
//            els_ip3derr   -> push_back( -999. );
//        }

	
    //Impact Parameters
	els_ip3d   -> push_back( el->dB(pat::Electron::PV3D) ); 
	els_ip3derr-> push_back( el->edB(pat::Electron::PV3D) ); 
	els_ip2d   -> push_back( el->dB(pat::Electron::PV2D) ); 
	els_ip2derr-> push_back( el->edB(pat::Electron::PV2D) ); 
	els_bs3d   -> push_back( el->dB(pat::Electron::BS3D) ); 
	els_bs3derr-> push_back( el->edB(pat::Electron::BS3D) ); 
	els_bs2d   -> push_back( el->dB(pat::Electron::BS2D) ); 
	els_bs2derr-> push_back( el->edB(pat::Electron::BS2D) ); 
	//els_ip3d      -> push_back( el->ip3d() ); // miniAOD


        /////////////////
        // Hit Pattern //
        /////////////////

        if( el_track->extra().isAvailable() ) {
            els_inner_position ->push_back(LorentzVector(el_track->innerPosition().x(), el_track->innerPosition().y() , el_track->innerPosition().z(), 0 ));
            els_outer_position ->push_back(LorentzVector(el_track->outerPosition().x(), el_track->outerPosition().y() , el_track->outerPosition().z(), 0 ));
        } else {
            els_inner_position->push_back(LorentzVector(-9999., -9999., -9999., -9999.));
            els_outer_position->push_back(LorentzVector(-9999., -9999., -9999., -9999.));
        }
    
	// Redesign according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/TrackingHitPatternRedesign
        const HitPattern& pattern = el_track->hitPattern();
        //const HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 
        //const HitPattern& p_outer = el_track->trackerExpectedHitsOuter();

        els_exp_innerlayers -> push_back(pattern.numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
        els_exp_outerlayers -> push_back(pattern.numberOfHits(reco::HitPattern::MISSING_OUTER_HITS));
        els_valid_pixelhits -> push_back(pattern.numberOfValidPixelHits());
        els_lost_pixelhits  -> push_back(pattern.numberOfLostPixelHits(reco::HitPattern::TRACK_HITS)); // Not sure about this. Could be MISSING_INNER_HITS instead.
	els_nlayers         -> push_back(pattern.trackerLayersWithMeasurement());
	els_nlayers3D       -> push_back(pattern.pixelLayersWithMeasurement() + pattern.numberOfValidStripLayersWithMonoAndStereo());
	els_nlayersLost     -> push_back(pattern.trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS));

        if( el_track->extra().isAvailable() ) {

            bool valid_hit      = false;
            uint32_t hit_pattern; 
            int i_layer       = 1;
            //int side = -1;
            bool pixel_hit   = false;
            bool strip_hit   = false;
            //int pixel_sizeX;
            //int pixel_sizeY;
            //float pixel_charge;
            //int det;
            //int layer;

            for( trackingRecHit_iterator ihit = el_track->recHitsBegin(); ihit != el_track->recHitsEnd(); ++ihit ) { 

                if(i_layer > 1) break;

                int k       = ihit-el_track->recHitsBegin();
                hit_pattern = pattern.getHitPattern(reco::HitPattern::TRACK_HITS, k);
                valid_hit   = pattern.validHitFilter(hit_pattern);
                pixel_hit   = pattern.pixelHitFilter(hit_pattern);
                strip_hit   = pattern.stripHitFilter(hit_pattern);
                //side        = (int)pattern.getSide(hit_pattern);
                //det         = (int)pattern.getSubStructure(hit_pattern);
                //layer       = (int)pattern.getLayer(hit_pattern);

                if(!valid_hit) continue;

                if(pixel_hit){
        
                    const SiPixelRecHit *pixel_hit_cast = dynamic_cast<const SiPixelRecHit*>(&(**ihit));
                    assert(pixel_hit_cast != 0);
                    //pixel_ClusterRef const& pixel_cluster = pixel_hit_cast->cluster();

                    //pixel_sizeX  = (int)pixel_cluster->sizeX(); 
                    //pixel_sizeY  = (int)pixel_cluster->sizeY(); 
                    //pixel_charge = (float)pixel_cluster->charge();
        
                    if( i_layer == 1 ) {
                    //    els_layer1_sizerphi -> push_back(pixel_sizeX);
                    //    els_layer1_sizerz   -> push_back(pixel_sizeY);
                    //    els_layer1_charge   -> push_back(pixel_charge);
                    //    els_layer1_det      -> push_back(det);
                    //    els_layer1_layer    -> push_back(layer);
                        i_layer++;
                    }

                } // end pixel hit
        
                else if (strip_hit){

                    //
                    const SiStripRecHit1D *strip_hit_cast   = dynamic_cast<const SiStripRecHit1D*>(&(**ihit));
                    const SiStripRecHit2D *strip2d_hit_cast = dynamic_cast<const SiStripRecHit2D*>(&(**ihit));
                    ClusterRef cluster;
                    if(strip_hit_cast == NULL){
                        cluster = strip2d_hit_cast->cluster();
                    }
                    else { 
                        cluster = strip_hit_cast->cluster();
                    }        

                    //
                    int cluster_size   = (int)cluster->amplitudes().size();
                    int cluster_charge = 0;
                    //int max_strip_i    = max_element(cluster->amplitudes().begin(),cluster->amplitudes().end())-cluster->amplitudes().begin();
                    //double cluster_weight_size = 0.0;

                    for( int istrip = 0; istrip < cluster_size; istrip++ ){
                        cluster_charge += (int)cluster->amplitudes().at(istrip);
                        //cluster_weight_size += (istrip-max_strip_i)*(istrip-max_strip_i)*(cluster->amplitudes().at(istrip));
                    }
                    //cluster_weight_size = sqrt(cluster_weight_size/cluster_charge);
        
                    if( i_layer == 1 ) {

                    //    //
                    //    els_layer1_charge -> push_back(cluster_charge);
                    //    els_layer1_det    -> push_back(det);
                    //    els_layer1_layer  -> push_back(layer);

                    //    //
                    //    if( side == 0 ) {
                    //        els_layer1_sizerphi -> push_back(cluster_size);
                    //        els_layer1_sizerz   -> push_back(0);
                    //    }
                    //    else {
                    //        els_layer1_sizerphi -> push_back(0);
                    //        els_layer1_sizerz   -> push_back(cluster_size);
                    //    }

                        i_layer++;

                    } // end layer = 1

                } // end strip hit

            } // end for loop

        } // end if extra 
        //else {
        //    els_layer1_sizerphi -> push_back(-9999);
        //    els_layer1_sizerz   -> push_back(-9999);
        //    els_layer1_charge   -> push_back(-9999);
        //    els_layer1_det      -> push_back(-9999);
        //    els_layer1_layer    -> push_back(-9999);
        //}
    

//        /////////////////
//        // Conversions //
//        /////////////////
//
//        ConversionFinder convFinder; //vector of conversion infos - all the candidate conversion tracks
//        vector<ConversionInfo> v_convInfos = convFinder.getConversionInfos(*(el->core()), tracks_h, gsftracks_h, evt_bField);
//    
//        vector<int>           v_tkidx;
//        vector<int>           v_gsftkidx;
//        vector<int>           v_delmisshits;
//        vector<int>           v_flag;
//        vector<float>         v_dist;
//        vector<float>         v_dcot;
//        vector<float>         v_rad;
//        vector<LorentzVector> v_pos_p4;
//
//        for(unsigned int i_conv = 0; i_conv < v_convInfos.size(); i_conv++) {
//      
//            //
//            math::XYZPoint convPoint  = v_convInfos.at(i_conv).pointOfConversion();
//            float          convPointX = isfinite(convPoint.x()) ? convPoint.x() : -9999.;
//            float          convPointY = isfinite(convPoint.y()) ? convPoint.y() : -9999.;
//            float          convPointZ = isfinite(convPoint.z()) ? convPoint.z() : -9999.;
//
//            //
//            v_dist        .push_back( isfinite(v_convInfos.at(i_conv).dist()) ? v_convInfos.at(i_conv).dist() : -9999.  );
//            v_dcot        .push_back( v_convInfos.at(i_conv).dcot()                                                     );
//            v_rad         .push_back( v_convInfos.at(i_conv).radiusOfConversion()                                       );
//            v_delmisshits .push_back( v_convInfos.at(i_conv).deltaMissingHits()                                         );
//            v_flag        .push_back( v_convInfos.at(i_conv).flag()                                                     );
//            v_pos_p4      .push_back( LorentzVector(convPointX, convPointY, convPointZ, 0)                              );
//
//            //
//            if( v_convInfos.at(i_conv).conversionPartnerCtfTk().isNonnull() ) {
//                v_tkidx.push_back(v_convInfos.at(i_conv).conversionPartnerCtfTk().key());
//            }
//            else {
//                v_tkidx.push_back(-9999);
//            }
//
//            //
//            if( v_convInfos.at(i_conv).conversionPartnerGsfTk().isNonnull() ) {
//                v_gsftkidx.push_back(v_convInfos.at(i_conv).conversionPartnerGsfTk().key());
//            }
//            else { 
//                v_gsftkidx.push_back(-9999);
//            }
//
//        } // end for loop
//
//
//        //
//        els_convs_dist->push_back(v_dist);
//        els_convs_dcot->push_back(v_dcot);
//        els_convs_radius->push_back(v_rad);
//        els_convs_tkidx->push_back(v_tkidx);
//        els_convs_gsftkidx->push_back(v_gsftkidx);
//        els_convs_delMissHits->push_back(v_delmisshits);
//        els_convs_flag->push_back(v_flag);
//
//        //
//        ConversionInfo convInfo   = convFinder.getConversionInfo( *el, tracks_h, gsftracks_h, evt_bField );
//        math::XYZPoint convPoint  = convInfo.pointOfConversion();
//        float          convPointX = isfinite(convPoint.x()) ? convPoint.x() : -9999.;
//        float          convPointY = isfinite(convPoint.y()) ? convPoint.y() : -9999.;
//        float          convPointZ = isfinite(convPoint.z()) ? convPoint.z() : -9999.;
//
//        //
//        els_conv_dist        -> push_back( isfinite(convInfo.dist()) ? convInfo.dist() : -9999. );
//        els_conv_dcot        -> push_back( convInfo.dcot()                                      );
//        els_conv_radius      -> push_back( convInfo.radiusOfConversion()                        );
//        els_conv_delMissHits -> push_back( convInfo.deltaMissingHits()                          );
//        els_conv_flag        -> push_back( convInfo.flag()                                      );
//        els_conv_pos_p4      -> push_back( LorentzVector(convPointX, convPointY, convPointZ, 0) );
//
//        //
//        if( convInfo.conversionPartnerCtfTk().isNonnull() ) {
//            els_conv_tkidx->push_back(convInfo.conversionPartnerCtfTk().key());
//        }
//        else {
//            els_conv_tkidx->push_back(-9999);
//        }
//
//        //
//        if( convInfo.conversionPartnerGsfTk().isNonnull() ) {
//            els_conv_gsftkidx->push_back(convInfo.conversionPartnerGsfTk().key());
//        }
//        else { 
//            els_conv_gsftkidx->push_back(-9999);
//        }


        //////////////////////////////
        // Flag For Vertex Fit Conversion Rejection //
        //////////////////////////////

//        Handle<ConversionCollection> convs_h;
//        iEvent.getByLabel(recoConversionInputTag_, convs_h);
//        els_conv_vtx_flag        -> push_back( ConversionTools::hasMatchedConversion(*el, convs_h, beamSpot ) );
        els_conv_vtx_flag        -> push_back( !el->passConversionVeto() ); // PAT variable: http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatAlgos/plugins/PATElectronProducer.cc#467


        ///////////////////
        // Added for 53x //
        ///////////////////

        els_passingMvaPreselection  ->push_back( el->passingMvaPreselection()   );
        els_passingPflowPreselection->push_back( el->passingPflowPreselection() );
        els_r9                      ->push_back( el->r9()                       );

        ///////////////////
        // Added for 7   //
        ///////////////////

        els_mass                    ->push_back( el->mass()                     );

	// Loop over PF candidates and find those associated by the map to the gedGsfElectron1
	vector<int> v_PFCand_idx;
	for( const edm::Ref<pat::PackedCandidateCollection> & ref : el->associatedPackedPFCandidates() )
	  v_PFCand_idx.push_back(ref.key());
	els_PFCand_idx->push_back(v_PFCand_idx);

	//////////////////////
	// genMatch miniAOD //
	//////////////////////
	
	LorentzVector mc_p4(0,0,0,0);	 
	const reco::GenParticle * gen = el->genParticle();
	if (gen != 0) {
	  mc_p4 = gen->p4();
	  els_mc_patMatch_id      ->push_back( gen->pdgId()  );
	  els_mc_patMatch_p4      ->push_back( mc_p4         );
	  els_mc_patMatch_dr      ->push_back( ROOT::Math::VectorUtil::DeltaR(gen->p4(), el->p4())  );
	}
	else {
	  els_mc_patMatch_id      ->push_back( -999   );
	  els_mc_patMatch_p4      ->push_back( mc_p4  );
	  els_mc_patMatch_dr      ->push_back( -999.  );
	}

	//////////////////////
	// mini-isolation   //
	//////////////////////

	float minichiso     = 0.;
	float mininhiso     = 0.;
	float miniemiso     = 0.;
	float minidbiso     = 0.;
	elMiniIso(el, true, 0.0, minichiso, mininhiso, miniemiso, minidbiso);
	els_miniIso_uncor   ->push_back( minichiso + mininhiso + miniemiso );
	els_miniIso_ch      ->push_back( minichiso );
	els_miniIso_nh      ->push_back( mininhiso );
	els_miniIso_em      ->push_back( miniemiso );
	els_miniIso_db      ->push_back( minidbiso );

	///////////////////////////
	// PFCluster isolation   //
	///////////////////////////

	els_ecalPFClusterIso ->push_back(   el->ecalPFClusterIso()  );
	els_hcalPFClusterIso ->push_back(   el->hcalPFClusterIso()  );



    } // end Loop on Electrons
  





    // Put the results into the event
    //
    iEvent.put(evt_nels, "evtnels");

    // Predefined ID descisions 
    //
    iEvent.put(els_class    , "elsclass"    );
    iEvent.put(els_category , "elscategory" );
  
    iEvent.put(passVetoId,   "passVetoId"   );
    iEvent.put(passLooseId,  "passLooseId"  );
    iEvent.put(passMediumId, "passMediumId" );
    iEvent.put(passTightId,  "passTightId"  );

    // Track parameters
    //
    iEvent.put(els_d0         , "elsd0"        );
    iEvent.put(els_z0         , "elsz0"        );
    iEvent.put(els_dxyPV      , "elsdxyPV"    );
    iEvent.put(els_dzPV       , "elsdzPV"     );
    iEvent.put(els_d0corr     , "elsd0corr"    );
    iEvent.put(els_d0corrPhi  , "elsd0corrPhi" );
    iEvent.put(els_d0phiCov   , "elsd0phiCov"  );
    iEvent.put(els_z0corr     , "elsz0corr"    );
    iEvent.put(els_chi2       , "elschi2"      );
    iEvent.put(els_ndof       , "elsndof"      );
    iEvent.put(els_d0Err      , "elsd0Err"     );
    iEvent.put(els_z0Err      , "elsz0Err"     );
    iEvent.put(els_ptErr      , "elsptErr"     );
    iEvent.put(els_ptErrGsf   , "elsptErrGsf"  );
    iEvent.put(els_etaErr     , "elsetaErr"    );
    iEvent.put(els_phiErr     , "elsphiErr"    );
    iEvent.put(els_ip3d       , "elsip3d"      );
    iEvent.put(els_ip3derr    , "elsip3derr"   );
    iEvent.put(els_ip2d       , "elsip2d"      );
    iEvent.put(els_ip2derr    , "elsip2derr"   );
    iEvent.put(els_bs3d       , "elsbs3d"      );
    iEvent.put(els_bs3derr    , "elsbs3derr"   );
    iEvent.put(els_bs2d       , "elsbs2d"      );
    iEvent.put(els_bs2derr    , "elsbs2derr"   );
  
    iEvent.put(els_validHits  , "elsvalidHits" );
    iEvent.put(els_lostHits   , "elslostHits"  );
    iEvent.put(els_charge     , "elscharge"    );
    iEvent.put(els_trk_charge , "elstrkcharge" );
    iEvent.put(els_sccharge   , "elssccharge"  );

    // Supercluster parameters
    //
    iEvent.put(els_nSeed       , "elsnSeed"       );
    iEvent.put(els_etaSC       , "elsetaSC"       );
    iEvent.put(els_phiSC       , "elsphiSC"       );
    iEvent.put(els_eSC         , "elseSC"         );
    iEvent.put(els_eSCRaw      , "elseSCRaw"      );
    iEvent.put(els_eSCPresh    , "elseSCPresh"    );
    iEvent.put(els_e1x5        , "else1x5"        );
//    iEvent.put(els_e3x3        , "else3x3"        );
    iEvent.put(els_e5x5        , "else5x5"        );
    iEvent.put(els_e2x5Max     , "else2x5Max"     );
//    iEvent.put(els_eMax        , "elseMax"        );
    iEvent.put(els_eSeed       , "elseSeed"       );
    iEvent.put(els_fiduciality , "elsfiduciality" );
    iEvent.put(els_type        , "elstype"        );
//    iEvent.put(els_scindex     , "elsscindex"     );
    iEvent.put(els_etaSCwidth  , "elsetaSCwidth"  );
    iEvent.put(els_phiSCwidth  , "elsphiSCwidth"  );

    // Corrections and uncertainties
    //
    iEvent.put(els_ecalEnergy         , "elsecalEnergy"         );
    iEvent.put(els_ecalEnergyError    , "elsecalEnergyError"    );
    iEvent.put(els_trackMomentumError , "elstrackMomentumError" );

    // Electron ID
    //
    //iEvent.put(els_sigmaPhiPhi        , "elssigmaPhiPhi"        );
    iEvent.put(els_sigmaIPhiIPhi      , "elssigmaIPhiIPhi"      );
//    iEvent.put(els_sigmaIEtaIPhi      , "elssigmaIEtaIPhi"      );
    iEvent.put(els_sigmaEtaEta        , "elssigmaEtaEta"        );
    iEvent.put(els_sigmaIEtaIEta      , "elssigmaIEtaIEta"      );
    //iEvent.put(els_sigmaIPhiIPhiSC    , "elssigmaIPhiIPhiSC"    );
    //iEvent.put(els_sigmaIEtaIEtaSC    , "elssigmaIEtaIEtaSC"    );
    iEvent.put(els_dPhiInPhiOut       , "elsdPhiInPhiOut"       );
    iEvent.put(els_hOverE             , "elshOverE"             );
    iEvent.put(els_hOverEBC           , "elshOverEBC"           );
    iEvent.put(els_hcalDepth1OverEcal , "elshcalDepth1OverEcal" );
    iEvent.put(els_hcalDepth2OverEcal , "elshcalDepth2OverEcal" );

    iEvent.put(els_eOverPIn                      , "elseOverPIn"                      );
    iEvent.put(els_eSeedOverPOut                 , "elseSeedOverPOut"                 );
    iEvent.put(els_eSeedOverPIn                  , "elseSeedOverPIn"                  );
    iEvent.put(els_eOverPOut                     , "elseOverPOut"                     );
    iEvent.put(els_fbrem                         , "elsfbrem"                         );
    //iEvent.put(els_lh                            , "elslh"                            );
    iEvent.put(els_mva                           , "elsmva"                           );
    iEvent.put(els_dEtaIn                        , "elsdEtaIn"                        );
    iEvent.put(els_dEtaOut                       , "elsdEtaOut"                       );
    iEvent.put(els_deltaEtaEleClusterTrackAtCalo , "elsdeltaEtaEleClusterTrackAtCalo" );
    iEvent.put(els_dPhiIn                        , "elsdPhiIn"                        );
    iEvent.put(els_dPhiOut                       , "elsdPhiOut"                       );
    iEvent.put(els_deltaPhiEleClusterTrackAtCalo , "elsdeltaPhiEleClusterTrackAtCalo" );
    iEvent.put(els_isGsfCtfScPixChargeConsistent , "elsisGsfCtfScPixChargeConsistent" );

    // Lorentz vectors
    //
    iEvent.put(els_p4     , "elsp4"    );
    iEvent.put(els_trk_p4 , "elstrkp4" );
    iEvent.put(els_p4In   , "elsp4In"  );
    iEvent.put(els_p4Out  , "elsp4Out" );

    // Vertex
    //
    iEvent.put(els_vertex_p4, "elsvertexp4");
    iEvent.put(els_trk_vertex_p4, "elstrkvertexp4");

    // Isolation
    //
    iEvent.put(els_tkIso                , "elstkIso"                );
    iEvent.put(els_ecalIso              , "elsecalIso"              );
    iEvent.put(els_hcalIso              , "elshcalIso"              );
    iEvent.put(els_hcalDepth1TowerSumEt , "elshcalDepth1TowerSumEt" );
    iEvent.put(els_hcalDepth2TowerSumEt , "elshcalDepth2TowerSumEt" );

    iEvent.put(els_tkIso04                , "elstkIso04"                );
    iEvent.put(els_ecalIso04              , "elsecalIso04"              );
    iEvent.put(els_hcalIso04              , "elshcalIso04"              );
    iEvent.put(els_hcalDepth1TowerSumEt04 , "elshcalDepth1TowerSumEt04" );
    iEvent.put(els_hcalDepth2TowerSumEt04 , "elshcalDepth2TowerSumEt04" );

//    iEvent.put(els_iso03_pf, "elsiso03pf" );
//    iEvent.put(els_iso04_pf, "elsiso04pf" );

//    iEvent.put(els_iso03_pf_ch      , "elsiso03pfch"      );
//    iEvent.put(els_iso03_pf_gamma05 , "elsiso03pfgamma05" );
//    iEvent.put(els_iso03_pf_nhad05  , "elsiso03pfnhad05"  );
//    iEvent.put(els_iso04_pf_ch      , "elsiso04pfch"      );
//    iEvent.put(els_iso04_pf_gamma05 , "elsiso04pfgamma05" );
//    iEvent.put(els_iso04_pf_nhad05  , "elsiso04pfnhad05"  );

    iEvent.put(els_pfChargedHadronIso , "elspfChargedHadronIso" );
    iEvent.put(els_pfNeutralHadronIso , "elspfNeutralHadronIso" );
    iEvent.put(els_pfPhotonIso        , "elspfPhotonIso"        );
    iEvent.put(els_pfPUIso            , "elspfPUIso"            );

    //Hit Pattern Information
    iEvent.put(els_inner_position  , "elsinnerposition"  );
    iEvent.put(els_outer_position  , "elsouterposition"  );
    iEvent.put(els_valid_pixelhits , "elsvalidpixelhits" );
    iEvent.put(els_lost_pixelhits  , "elslostpixelhits"  );
    iEvent.put(els_nlayers         , "elsnlayers"        );
    iEvent.put(els_nlayers3D       , "elsnlayers3D"      );
    iEvent.put(els_nlayersLost     , "elsnlayersLost"    );
    //iEvent.put(els_layer1_layer    , "elslayer1layer"    );
    //iEvent.put(els_layer1_sizerphi , "elslayer1sizerphi" );
    //iEvent.put(els_layer1_sizerz   , "elslayer1sizerz"   );
    //iEvent.put(els_layer1_charge   , "elslayer1charge"   );
    //iEvent.put(els_layer1_det      , "elslayer1det"      );
    iEvent.put(els_exp_innerlayers , "elsexpinnerlayers" );
    iEvent.put(els_exp_outerlayers , "elsexpouterlayers" );

    //CTF track info
    //
    iEvent.put(els_trkdr           , "elstrkdr"         );
    iEvent.put(els_trkshFrac       , "elstrkshFrac"     );
    iEvent.put(els_ckf_chi2        ,"elsckfchi2"        );
    iEvent.put(els_ckf_ndof        ,"elsckfndof"        );
    iEvent.put(els_ckf_laywithmeas ,"elsckflaywithmeas" );
    iEvent.put(els_ckf_charge      ,"elsckfcharge"      );

    //conversion
    //iEvent.put(els_convs_dist        , "elsconvsdist"        );
    //iEvent.put(els_convs_dcot        , "elsconvsdcot"        );
    //iEvent.put(els_convs_radius      , "elsconvsradius"      );
    //iEvent.put(els_convs_tkidx       , "elsconvstkidx"       );
    //iEvent.put(els_convs_gsftkidx    , "elsconvsgsftkidx"    );
    //iEvent.put(els_convs_delMissHits , "elsconvsdelMissHits" );
    //iEvent.put(els_convs_flag        , "elsconvsflag"        );


    //iEvent.put(els_conv_dist        , "elsconvdist"        );
    //iEvent.put(els_conv_dcot        , "elsconvdcot"        );
    //iEvent.put(els_conv_radius      , "elsconvradius"      );
    //iEvent.put(els_conv_pos_p4      , "elsconvposp4"       );
    //iEvent.put(els_conv_tkidx       , "elsconvtkidx"       );
    //iEvent.put(els_conv_gsftkidx    , "elsconvgsftkidx"    );
    //iEvent.put(els_conv_delMissHits , "elsconvdelMissHits" );
    //iEvent.put(els_conv_flag        , "elsconvflag"        );

    //iEvent.put(els_conv_old_dist        , "elsconvolddist"        );
    //iEvent.put(els_conv_old_dcot        , "elsconvolddcot"        );
    //iEvent.put(els_conv_old_radius      , "elsconvoldradius"      );
    //iEvent.put(els_conv_old_tkidx       , "elsconvoldtkidx"       );
    //iEvent.put(els_conv_old_gsftkidx    , "elsconvoldgsftkidx"    );
    //iEvent.put(els_conv_old_delMissHits , "elsconvolddelMissHits" );
    //iEvent.put(els_conv_old_flag        , "elsconvoldflag"        );

    iEvent.put(els_conv_vtx_flag        , "elsconvvtxflag"        );

    ///////////////////
    // Added for 53x //
    ///////////////////

    iEvent.put( els_passingMvaPreselection   , "elspassingMvaPreselection"   );
    iEvent.put( els_passingPflowPreselection , "elspassingPflowPreselection" );
    iEvent.put( els_r9                       , "elsr9"                       );

    ///////////////////
    // Added for 7   //
    ///////////////////

    iEvent.put(els_PFCand_idx    , "elspfcandidx"    );
    iEvent.put(els_mass          , "elsmass"         );

    /////////////////////////
    // Added for miniAOD   //
    /////////////////////////

    // genParticle matching from miniAOD
  iEvent.put( els_mc_patMatch_id        , "elsmcpatMatchid"        	);
  iEvent.put( els_mc_patMatch_p4        , "elsmcpatMatchp4"         );
  iEvent.put( els_mc_patMatch_dr        , "elsmcpatMatchdr"         );
  iEvent.put(els_sigmaIPhiIPhi_full5x5  , "elssigmaIPhiIPhifull5x5" );
  iEvent.put(els_sigmaEtaEta_full5x5    , "elssigmaEtaEtafull5x5"   );
  iEvent.put(els_sigmaIEtaIEta_full5x5  , "elssigmaIEtaIEtafull5x5" );
  iEvent.put(els_r9_full5x5             , "elsr9full5x5"            );
  iEvent.put(els_e1x5_full5x5           , "else1x5full5x5"          );
  iEvent.put(els_e5x5_full5x5           , "else5x5full5x5"          );
  iEvent.put(els_e2x5Max_full5x5        , "else2x5Maxfull5x5"       ); 

  iEvent.put(els_miniIso_uncor       , "elsminiIsouncor"    );
  iEvent.put(els_miniIso_ch       , "elsminiIsoch"    );
  iEvent.put(els_miniIso_nh       , "elsminiIsonh"    );
  iEvent.put(els_miniIso_em       , "elsminiIsoem"    );
  iEvent.put(els_miniIso_db       , "elsminiIsodb"    );

  iEvent.put(els_ecalPFClusterIso       , "elsecalPFClusterIso"    );
  iEvent.put(els_hcalPFClusterIso       , "elshcalPFClusterIso"    );
  
}

//----------------------------------------------------------------------------
// Electron Id classification function (a flag for the Sani type class)
//----------------------------------------------------------------------------
int ElectronMaker::classify(const RefToBase<pat::Electron> &electron) {

    double eOverP = electron->eSuperClusterOverP();
    double fbrem = electron->fbrem();
  
    int cat;
    if((electron->isEB() && fbrem<0.06) || (electron->isEE() && fbrem<0.1)) 
        cat=1;
    else if (eOverP < 1.2 && eOverP > 0.8) 
        cat=0;
    else 
        cat=2;
  
    return cat;

}

//little labour saving function to get the reference to the ValueMap
template<typename T> const ValueMap<T>& ElectronMaker::getValueMap(const Event& iEvent, InputTag& inputTag){
    Handle<ValueMap<T> > handle;
    iEvent.getByLabel(inputTag,handle);
    return *(handle.product());
}

double ElectronMaker::electronIsoValuePF(const GsfElectron& el, const Vertex& vtx, float coner, float minptn, float dzcut, float footprintdr, float gammastripveto, float elestripveto, int filterId){

    float pfciso = 0.;
    float pfniso = 0.;
    float pffootprint = 0.;
    float pfjurveto = 0.;
    float pfjurvetoq = 0.;

    //TrackRef siTrack     = el.closestCtfTrackRef();
    TrackRef siTrack     = el.closestTrack();
    GsfTrackRef gsfTrack = el.gsfTrack();

    if (gsfTrack.isNull() && siTrack.isNull()) return -9999.;

    float eldz = gsfTrack.isNonnull() ? gsfTrack->dz(vtx.position()) : siTrack->dz(vtx.position());
    float eleta = el.eta();

    for (PFCandidateCollection::const_iterator pf=pfCand_h->begin(); pf<pfCand_h->end(); ++pf){

        float pfeta = pf->eta();    
        float dR = deltaR(pfeta, pf->phi(), eleta, el.phi());
        if (dR>coner) continue;

        float deta = fabs(pfeta - eleta);
        int pfid = abs(pf->pdgId());
        float pfpt = pf->pt();

        if (filterId!=0 && filterId!=pfid) continue;

        if (pf->charge()==0) {
            //neutrals
            if (pfpt>minptn) {
                pfniso+=pfpt;
                if (dR<footprintdr && pfid==130) pffootprint+=pfpt;
                if (deta<gammastripveto && pfid==22)  pfjurveto+=pfpt;
            }
        } else {
            //charged  
            //avoid double counting of electron itself
            //if either the gsf or the ctf track are shared with the candidate, skip it
            const TrackRef pfTrack  = pf->trackRef();
            if (siTrack.isNonnull()  && pfTrack.isNonnull() && siTrack.key()==pfTrack.key()) continue;
            //below pfid==1 is commented out: in some cases the pfCand has a gsf even if it is not an electron... this is to improve the sync with MIT
            if (/*pfid==11 &&*/ pf->gsfTrackRef().isNonnull()) {
                if (gsfTrack.isNonnull() && gsfTrack.key()==pf->gsfTrackRef().key()) continue;
            } 
            //check electrons with gsf track
            if (pfid==11 && pf->gsfTrackRef().isNonnull()) {
                if(fabs(pf->gsfTrackRef()->dz(vtx.position()) - eldz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                    if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
                }
                continue;//and avoid double counting
            }
            //then check anything that has a ctf track
            if (pfTrack.isNonnull()) {//charged (with a ctf track)
                if(fabs( pfTrack->dz(vtx.position()) - eldz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                    if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
                }
            }
        } 
    }
    return pfciso+pfniso-pffootprint-pfjurveto-pfjurvetoq;
}

void ElectronMaker::PFIsolation2012(const reco::GsfElectron& el, const reco::VertexCollection* vertexCollection,
        const int vertexIndex, const float &R, float &pfiso_ch, float &pfiso_em, float &pfiso_nh)
{

    // isolation sums
    pfiso_ch = 0.0;
    pfiso_em = 0.0;
    pfiso_nh = 0.0;

    // loop on pfcandidates
    reco::PFCandidateCollection::const_iterator pf = pfCand_h->begin();
    for (pf = pfCand_h->begin(); pf != pfCand_h->end(); ++pf) {

        // skip electrons and muons
        if (pf->particleId() == reco::PFCandidate::e)     continue;
        if (pf->particleId() == reco::PFCandidate::mu)    continue;

        // deltaR between electron and cadidate
        const float dR = deltaR(pf->eta(), pf->phi(), el.eta(), el.phi());
        if (dR > R)                             continue;

        // charged hadrons closest vertex
        // should be the primary vertex
        if (pf->particleId() == reco::PFCandidate::h) {
            int pfVertexIndex = pfPileUpAlgo_->chargedHadronVertex(*vertexCollection, *pf);
            if (pfVertexIndex != vertexIndex) continue;
        }

        // endcap region
        if (!el.isEB()) {
            if (pf->particleId() == reco::PFCandidate::h      && dR <= 0.015)   continue;
            if (pf->particleId() == reco::PFCandidate::gamma  && dR <= 0.08)    continue;
        }

        // add to isolation sum
        if (pf->particleId() == reco::PFCandidate::h)       pfiso_ch += pf->pt();
        if (pf->particleId() == reco::PFCandidate::gamma)   pfiso_em += pf->pt();
        if (pf->particleId() == reco::PFCandidate::h0)      pfiso_nh += pf->pt();

    }

}

void ElectronMaker::elIsoCustomCone(edm::View<pat::Electron>::const_iterator& el, float dr, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float & dbiso){
  chiso     = 0.;
  nhiso     = 0.;
  emiso     = 0.;
  dbiso     = 0.;
  float deadcone_ch = 0.;
  float deadcone_pu = 0.;
  float deadcone_ph = 0.;
  // veto cones only in the endcap for electrons
  if (useVetoCones && fabs(el->eta()) > 1.479) { 
    deadcone_ch = 0.015;
    deadcone_pu = 0.015;
    deadcone_ph = 0.08;
  }
  for( pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++ ) {
    float thisDR = fabs(ROOT::Math::VectorUtil::DeltaR(pf_it->p4(),el->p4()));
    if ( thisDR>dr ) continue;  
    float pt = pf_it->p4().pt();
    float id = pf_it->pdgId();
    if ( fabs(id)==211 ) {
      if (pf_it->fromPV() > 1 && (!useVetoCones || thisDR > deadcone_ch) ) chiso+=pt;
      else if ((pf_it->fromPV() <= 1) && (pt > ptthresh) && (!useVetoCones || thisDR > deadcone_pu)) dbiso+=pt;
    }
    if ( fabs(id)==130 && (pt > ptthresh) ) nhiso+=pt;
    if ( fabs(id)==22 && (pt > ptthresh) && (!useVetoCones || thisDR > deadcone_ph) ) emiso+=pt;
  }
  //if (useDBcor) correction = 0.5 * deltaBeta;
  //else if (useEAcor) correction = evt_fixgrid_all_rho() * elEA03(elIdx) * (dr/0.3) * (dr/0.3);
  //float absiso = chiso + std::max(float(0.0), nhiso + emiso - correction);
  return;
}

void ElectronMaker::elMiniIso(edm::View<pat::Electron>::const_iterator& el, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float &dbiso){

  float pt = el->p4().pt();
  float dr = 0.2;
  if (pt>50) dr = 10./pt;
  if (pt>200) dr = 0.05;
  elIsoCustomCone(el,dr,useVetoCones,ptthresh, chiso, nhiso, emiso, dbiso);
  return;
}


//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMaker);
