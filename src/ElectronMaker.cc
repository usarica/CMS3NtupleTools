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
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronMaker.cc,v 1.65 2011/04/28 00:50:29 dbarge Exp $
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
#include "FWCore/Framework/interface/ESHandle.h"

#include "CMS2/NtupleMaker/interface/ElectronMaker.h"
#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"
#include "CMS2/NtupleMaker/interface/EgammaFiduciality.h"

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
#include "Math/VectorUtil.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

//conversion
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// class decleration
//

//
// constructors and destructor
//
ElectronMaker::ElectronMaker(const edm::ParameterSet& iConfig) {

  //get setup parameters
  electronsInputTag_             = iConfig.getParameter<edm::InputTag>("electronsInputTag"                  );
  beamSpotInputTag_              = iConfig.getParameter<edm::InputTag>("beamSpotInputTag"                   );
  trksInputTag_                  = iConfig.getParameter<edm::InputTag>("trksInputTag"                       );

  gsftracksInputTag_             = iConfig.getParameter<edm::InputTag>("gsftracksInputTag"                  );
  cms2scsseeddetidInputTag_      = iConfig.getParameter<edm::InputTag>("cms2scsseeddetidInputTag");
  
  minAbsDist_                    = iConfig.getParameter<double>("minAbsDist"                         );
  minAbsDcot_                    = iConfig.getParameter<double>("minAbsDcot"                         );
  minSharedFractionOfHits_       = iConfig.getParameter<double>("minSharedFractionOfHits"            );
  
  mtsTransform_ = 0;
  clusterTools_ = 0;
  
  aliasprefix_            = iConfig.getUntrackedParameter<std::string>("aliasPrefix");


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
  produces<vector<int> >       ("elsscindex"                 ).setBranchAlias("els_scindex"                );

  // Corrections and uncertainties
  //
  produces<vector<float> >     ("elsecalEnergy"              ).setBranchAlias("els_ecalEnergy"             );
  produces<vector<float> >     ("elsecalEnergyError"         ).setBranchAlias("els_ecalEnergyError"        );
  produces<vector<float> >     ("elstrackMomentumError"      ).setBranchAlias("els_trackMomentumError"     );
  produces<vector<float> >     ("elselectronMomentumError"   ).setBranchAlias("els_electronMomentumError"  );

  // ID variables
  //
  produces<vector<float> >     ("elsmva"                  ).setBranchAlias("els_mva"                 );
  produces<vector<float> >     ("elsdEtaIn"                  ).setBranchAlias("els_dEtaIn"                 );
  produces<vector<float> >     ("elsdEtaOut"                 ).setBranchAlias("els_dEtaOut"                );
  produces<vector<float> >     ("elsdeltaEtaEleClusterTrackAtCalo").setBranchAlias("els_deltaEtaEleClusterTrackAtCalo");
  produces<vector<float> >     ("elsdPhiIn"                  ).setBranchAlias("els_dPhiIn"                 );
  produces<vector<float> >     ("elsdPhiOut"                 ).setBranchAlias("els_dPhiOut"                );
  produces<vector<float> >     ("elsdeltaPhiEleClusterTrackAtCalo").setBranchAlias("els_deltaPhiEleClusterTrackAtCalo");
  produces<vector<float> >     ("elsdPhiInPhiOut"            ).setBranchAlias("els_dPhiInPhiOut"           );
  produces<vector<float> >     ("elsfbrem"                   ).setBranchAlias("els_fbrem"                  );
  produces<vector<float> >     ("elseSeed"                   ).setBranchAlias("els_eSeed"                  );
  produces<vector<float> >     ("elseOverPIn"                ).setBranchAlias("els_eOverPIn"               );
  produces<vector<float> >     ("elseSeedOverPOut"           ).setBranchAlias("els_eSeedOverPOut"          );
  produces<vector<float> >     ("elseSeedOverPIn"            ).setBranchAlias("els_eSeedOverPIn"           );
  produces<vector<float> >     ("elseOverPOut"               ).setBranchAlias("els_eOverPOut"              );

  produces<vector<float> >     ("elshOverE"                  ).setBranchAlias("els_hOverE"                 );
  produces<vector<float> >     ("elshcalDepth1OverEcal"      ).setBranchAlias("els_hcalDepth1OverEcal"     );
  produces<vector<float> >     ("elshcalDepth2OverEcal"      ).setBranchAlias("els_hcalDepth2OverEcal"     );

  produces<vector<float> >     ("elssigmaPhiPhi"             ).setBranchAlias("els_sigmaPhiPhi"            );
  produces<vector<float> >     ("elssigmaIPhiIPhi"           ).setBranchAlias("els_sigmaIPhiIPhi"          );
  produces<vector<float> >     ("elssigmaEtaEta"             ).setBranchAlias("els_sigmaEtaEta"            );
  produces<vector<float> >     ("elssigmaIEtaIEta"           ).setBranchAlias("els_sigmaIEtaIEta"          );
  produces<vector<float> >     ("elssigmaIPhiIPhiSC"         ).setBranchAlias("els_sigmaIPhiIPhiSC"        );
  produces<vector<float> >     ("elssigmaIEtaIEtaSC"         ).setBranchAlias("els_sigmaIEtaIEtaSC"        );
  produces<vector<float> >     ("else2x5Max"                 ).setBranchAlias("els_e2x5Max"                );
  produces<vector<float> >     ("else1x5"                    ).setBranchAlias("els_e1x5"                   );
  produces<vector<float> >     ("else5x5"                    ).setBranchAlias("els_e5x5"                   );
  produces<vector<float> >     ("else3x3"                    ).setBranchAlias("els_e3x3"                   );
  produces<vector<float> >     ("elseMax"                    ).setBranchAlias("els_eMax"                   );

  // predefined ID decisions
  // http://cmslxr.fnal.gov/lxr/source/DataFormats/EgammaCandidates/interface/GsfElectron.h
  produces<vector<int> >       ("elsclass"                   ).setBranchAlias("els_class"                  );
  // for the ID definitions, see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideElectronID
  // the decisions should be the SAME as the els_pat_*id branches made by PATElectronMaker
  produces<vector<int> >       ("elscategory"                ).setBranchAlias("els_category"               );

  // isolation variables
  //
  produces<vector<float> >     ("elstkIso"                   ).setBranchAlias("els_tkIso"                  );
  produces<vector<float> >     ("elsecalIso"                 ).setBranchAlias("els_ecalIso"                );
  produces<vector<float> >     ("elshcalIso"                 ).setBranchAlias("els_hcalIso"                );
  produces<vector<float> >     ("elshcalDepth1TowerSumEt"    ).setBranchAlias("els_hcalDepth1TowerSumEt"   );
  produces<vector<float> >     ("elshcalDepth2TowerSumEt"    ).setBranchAlias("els_hcalDepth2TowerSumEt"   );

  produces<vector<float> >     ("elstkIso04"                 ).setBranchAlias("els_tkIso04"                );
  produces<vector<float> >     ("elsecalIso04"               ).setBranchAlias("els_ecalIso04"              );
  produces<vector<float> >     ("elshcalIso04"               ).setBranchAlias("els_hcalIso04"              );
  produces<vector<float> >     ("elshcalDepth1TowerSumEt04"  ).setBranchAlias("els_hcalDepth1TowerSumEt04" );
  produces<vector<float> >     ("elshcalDepth2TowerSumEt04"  ).setBranchAlias("els_hcalDepth2TowerSumEt04" );

  // track variables
  //
  produces<vector<int> >       ("elscharge"                  ).setBranchAlias("els_charge"                 ); //candidate charge
  produces<vector<int> >       ("elssccharge"                ).setBranchAlias("els_sccharge"               );
  produces<vector<int> >       ("elstrkcharge"               ).setBranchAlias("els_trk_charge"             );
  produces<vector<float> >     ("elschi2"                    ).setBranchAlias("els_chi2"                   );
  produces<vector<float> >     ("elsndof"                    ).setBranchAlias("els_ndof"                   );
  produces<vector<int> >       ("elsvalidHits"               ).setBranchAlias("els_validHits"              ); //number of used hits in fit
  produces<vector<int> >       ("elslostHits"                ).setBranchAlias("els_lostHits"               ); //number of lost hits in fit
  produces<vector<float> >     ("elsd0"                      ).setBranchAlias("els_d0"                     );
  produces<vector<float> >     ("elsz0"                      ).setBranchAlias("els_z0"                     );
  produces<vector<float> >     ("elsd0Err"                   ).setBranchAlias("els_d0Err"                  );
  produces<vector<float> >     ("elsz0Err"                   ).setBranchAlias("els_z0Err"                  );
  produces<vector<float> >     ("elsd0corr"                  ).setBranchAlias("els_d0corr"                 );
  produces<vector<float> >     ("elsz0corr"                  ).setBranchAlias("els_z0corr"                 );
  produces<vector<float> >     ("elsptErr"                   ).setBranchAlias("els_ptErr"                  );
  produces<vector<float> >     ("elsetaErr"                  ).setBranchAlias("els_etaErr"                 );
  produces<vector<float> >     ("elsphiErr"                  ).setBranchAlias("els_phiErr"                 );
  produces<vector<int> >                ("elsgsftrkidx"                         ).setBranchAlias("els_gsftrkidx"                    );
  
  // LorentzVectors
  //
  produces<vector<LorentzVector> >  ("elsp4"                 ).setBranchAlias("els_p4"                     );
  produces<vector<LorentzVector> >  ("elstrkp4"              ).setBranchAlias("els_trk_p4"                 );
  produces<vector<LorentzVector> >  ("elsp4In"               ).setBranchAlias("els_p4In"                   );
  produces<vector<LorentzVector> >  ("elsp4Out"              ).setBranchAlias("els_p4Out"                  );

  // Vertex
  //
  produces<vector<LorentzVector> >  ("elsvertexp4"           ).setBranchAlias("els_vertex_p4"              );

  //Hit Pattern information

  produces<vector<LorentzVector> >  ("elsinnerposition"           ).setBranchAlias("els_inner_position"         );
  produces<vector<LorentzVector> >  ("elsouterposition"           ).setBranchAlias("els_outer_position"         );
  produces<vector<int> >            ("elsvalidpixelhits"          ).setBranchAlias("els_valid_pixelhits"        );
  produces<vector<int> >            ("elslostpixelhits"           ).setBranchAlias("els_lost_pixelhits"         );
  produces<vector<int> >            ("elslayer1sizerphi"          ).setBranchAlias("els_layer1_sizerphi"        ); 
  produces<vector<int> >            ("elslayer1sizerz"            ).setBranchAlias("els_layer1_sizerz"          ); 
  produces<vector<float> >          ("elslayer1charge"            ).setBranchAlias("els_layer1_charge"          ); 
  produces<vector<int> >            ("elslayer1det"               ).setBranchAlias("els_layer1_det"             );
  produces<vector<int> >            ("elslayer1layer"             ).setBranchAlias("els_layer1_layer"           ); 
  produces<vector<int> >            ("elsexpinnerlayers"          ).setBranchAlias("els_exp_innerlayers"        );
  produces<vector<int> >            ("elsexpouterlayers"          ).setBranchAlias("els_exp_outerlayers"        );   
	
  //CTF track matching stuff
  produces<vector<int>    >    ("elstrkidx"                  ).setBranchAlias("els_trkidx"                 );// track index matched to electron
  produces<vector<float>  >    ("elstrkshFrac"               ).setBranchAlias("els_trkshFrac"              );
  produces<vector<float>  >    ("elstrkdr"                   ).setBranchAlias("els_trkdr"                  );


  //These are vectors of vectors, holding the candidate conversion partners. 
  produces<vector<vector<float> >   >       ("elsconvsdist"          ).setBranchAlias("els_convs_dist"            );
  produces<vector<vector<float> >   >       ("elsconvsdcot"          ).setBranchAlias("els_convs_dcot"            );
  produces<vector<vector<float> >   >       ("elsconvsradius"        ).setBranchAlias("els_convs_radius"          );  //signed radius of conversion
  produces<vector<vector<LorentzVector> > > ("elsconvsposp4"         ).setBranchAlias("els_convs_pos_p4"          );  //position of conversion
  produces<vector<vector<int>   >   >       ("elsconvstkidx"         ).setBranchAlias("els_convs_tkidx"           );  //index of partner track
  produces<vector<vector<int>   >   >       ("elsconvsgsftkidx"      ).setBranchAlias("els_convs_gsftkidx"        );  //index of the GSF partner track, if thats where we find it
  produces<vector<vector<int>   >   >       ("elsconvsdelMissHits"   ).setBranchAlias("els_convs_delMissHits"     );  //Delta Missing Hits between the electron and partner track
  produces<vector<vector<int>   >   >       ("elsconvsflag"          ).setBranchAlias("els_convs_flag"            );


  //conversion stuff - This is the "best" conversion partner, defined as the one that 
  //has the minimum sqrt(dist*dist + dcot*dcot)
  produces<vector<float>    >       ("elsconvdist"          ).setBranchAlias("els_conv_dist"            );
  produces<vector<float>    >       ("elsconvdcot"          ).setBranchAlias("els_conv_dcot"            );
  produces<vector<float>    >       ("elsconvradius"        ).setBranchAlias("els_conv_radius"          );  //signed radius of conversion
  produces<vector<LorentzVector> >  ("elsconvposp4"         ).setBranchAlias("els_conv_pos_p4"          );  //position of conversion
  produces<vector<int>      >       ("elsconvtkidx"         ).setBranchAlias("els_conv_tkidx"           );  //index of partner track
  produces<vector<int>      >       ("elsconvgsftkidx"      ).setBranchAlias("els_conv_gsftkidx"        );  //index of the GSF partner track, if thats where we find it
  produces<vector<int>      >       ("elsconvdelMissHits"   ).setBranchAlias("els_conv_delMissHits"     );  //Delta Missing Hits between the electron and partner track
  produces<vector<int>      >       ("elsconvflag"          ).setBranchAlias("els_conv_flag"            );
  


  produces<vector<float>    >       ("elsconvolddist"       ).setBranchAlias("els_conv_old_dist"        );
  produces<vector<float>    >       ("elsconvolddcot"       ).setBranchAlias("els_conv_old_dcot"        );
  produces<vector<float>    >       ("elsconvoldradius"     ).setBranchAlias("els_conv_old_radius"      );  //signed radius of conversion
  produces<vector<int>      >       ("elsconvoldtkidx"      ).setBranchAlias("els_conv_old_tkidx"       );  //index of partner track
  produces<vector<int>      >       ("elsconvoldgsftkidx"   ).setBranchAlias("els_conv_old_gsftkidx"    );  //index of the GSF partner track, if thats where we find it
  produces<vector<int>      >       ("elsconvolddelMissHits").setBranchAlias("els_conv_old_delMissHits" );  //Delta Missing Hits between the electron and partner track
  produces<vector<int>      >       ("elsconvoldflag"       ).setBranchAlias("els_conv_old_flag"        );

  
}

ElectronMaker::~ElectronMaker()
{
  if (clusterTools_) delete clusterTools_;
  if (mtsTransform_) delete mtsTransform_;
}
/*
void  ElectronMaker::beginJob() {
  
  edm::ESHandle<TrackerGeometry>              trackerGeometryHandle;
  edm::ESHandle<MagneticField>                magFieldHandle;
  
  es.get<TrackerDigiGeometryRecord>().get(trackerGeometryHandle);
  es.get<IdealMagneticFieldRecord>().get(magFieldHandle);
  mtsTransform_ = new MultiTrajectoryStateTransform(trackerGeometryHandle.product(), magFieldHandle.product());
   

}
*/

void  ElectronMaker::beginRun(edm::Run&, const edm::EventSetup& es) {
  
  edm::ESHandle<TrackerGeometry>              trackerGeometryHandle;
  edm::ESHandle<MagneticField>                magFieldHandle;
  
  es.get<TrackerDigiGeometryRecord>().get(trackerGeometryHandle);
  es.get<IdealMagneticFieldRecord>().get(magFieldHandle);
  mtsTransform_ = new MultiTrajectoryStateTransform(trackerGeometryHandle.product(), magFieldHandle.product());
   

}

void ElectronMaker::beginJob() {
}

void ElectronMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void ElectronMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Define vectors to be filled

  auto_ptr<unsigned int>	   evt_nels				(new unsigned int        ) ;

  // ECAL related (superCluster) variables
  auto_ptr<vector<int> >		els_nSeed				(new vector<int>		) ;
  auto_ptr<vector<float> >		els_etaSC				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_phiSC				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_eSC					(new vector<float>		) ;
  auto_ptr<vector<float> >		els_eSCRaw				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_eSCPresh				(new vector<float>		) ;
  auto_ptr<vector<int> >		els_fiduciality				(new vector<int>		) ;
  auto_ptr<vector<int> >		els_type				(new vector<int>		) ;
  auto_ptr<vector<int> >		els_scindex				(new vector<int>		) ;	 

  // uncertainties and corrections
  // somewhat complicated: see 
  // http://cms-service-sdtweb.web.cern.ch/cms-service-sdtweb/doxygen/CMSSW_3_1_2/doc/html/d5/d4b/GsfElectron_8h-source.html
  // note that if ecalEnergy == eSC depends on if further ecal corrections have been applied to the electron
  // after its construction
  auto_ptr<vector<float> >		els_ecalEnergy				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_ecalEnergyError			(new vector<float>		) ;
  auto_ptr<vector<float> >		els_trackMomentumError			(new vector<float>		) ;
  auto_ptr<vector<float> >		els_electronMomentumError		(new vector<float>		) ;
  
  // ID variables
  //
  auto_ptr<vector<float> >		els_mva					(new vector<float>		) ;
  auto_ptr<vector<float> >		els_dEtaIn				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_dEtaOut				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_dPhiIn				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_dPhiOut				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_dPhiInPhiOut			(new vector<float>		) ;
  auto_ptr<vector<float> >		els_fbrem				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_eSeed				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_eOverPIn				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_eSeedOverPOut			(new vector<float>		) ;
  auto_ptr<vector<float> >		els_eSeedOverPIn			(new vector<float>		) ;
  auto_ptr<vector<float> >		els_eOverPOut				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_deltaEtaEleClusterTrackAtCalo	(new vector<float>		) ;
  auto_ptr<vector<float> >		els_deltaPhiEleClusterTrackAtCalo	(new vector<float>		) ;

  auto_ptr<vector<float> >		els_hOverE				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_hcalDepth1OverEcal			(new vector<float>		) ;
  auto_ptr<vector<float> >		els_hcalDepth2OverEcal			(new vector<float>		) ;

  auto_ptr<vector<float> >		els_sigmaPhiPhi				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_sigmaIPhiIPhi			(new vector<float>		) ;
  auto_ptr<vector<float> >		els_sigmaEtaEta				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_sigmaIEtaIEta			(new vector<float>		) ;
  auto_ptr<vector<float> >		els_sigmaIPhiIPhiSC			(new vector<float>		) ;
  auto_ptr<vector<float> >		els_sigmaIEtaIEtaSC			(new vector<float>		) ;

  auto_ptr<vector<float> >		els_e2x5Max				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_e1x5				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_e5x5				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_e3x3				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_eMax				(new vector<float>		) ;

  // predefined ID decisions
  //
  auto_ptr<vector<int> >		els_class				(new vector<int>		) ;
  auto_ptr<vector<int> >		els_category				(new vector<int>		) ;


  // isolation variables
  //
  auto_ptr<vector<float> >		els_tkIso				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_ecalIso				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_hcalIso				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_hcalDepth1TowerSumEt		(new vector<float>		) ;
  auto_ptr<vector<float> >		els_hcalDepth2TowerSumEt		(new vector<float>		) ;

  auto_ptr<vector<float> >		els_tkIso04				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_ecalIso04				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_hcalIso04				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_hcalDepth1TowerSumEt04		(new vector<float>		) ;
  auto_ptr<vector<float> >		els_hcalDepth2TowerSumEt04		(new vector<float>		) ;

  // track variables
  //
  auto_ptr<vector<int> >		els_charge				(new vector<int>		) ;
  auto_ptr<vector<int> >		els_trk_charge				(new vector<int>		) ;
  auto_ptr<vector<int> >		els_sccharge				(new vector<int>		) ;
  auto_ptr<vector<float> >		els_chi2				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_ndof				(new vector<float>		) ;
  auto_ptr<vector<int> >		els_validHits				(new vector<int>		) ;
  auto_ptr<vector<int> >		els_lostHits				(new vector<int>		) ;
  auto_ptr<vector<float> >		els_d0					(new vector<float>		) ;
  auto_ptr<vector<float> >		els_z0					(new vector<float>		) ;
  auto_ptr<vector<float> >		els_d0Err				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_z0Err				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_d0corr				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_z0corr				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_ptErr				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_etaErr				(new vector<float>		) ;
  auto_ptr<vector<float> >		els_phiErr				(new vector<float>		) ;
  auto_ptr<vector<int>   >		els_gsftrkidx				(new vector<int>                ) ;

  
  // LorentzVectors
  //
  auto_ptr<vector<LorentzVector> >	els_p4					(new vector<LorentzVector>	) ;
  auto_ptr<vector<LorentzVector> >	els_trk_p4				(new vector<LorentzVector>	) ;
  auto_ptr<vector<LorentzVector> >	els_p4In				(new vector<LorentzVector>	) ;
  auto_ptr<vector<LorentzVector> >	els_p4Out				(new vector<LorentzVector>	) ;

  // Vertex
  //
  auto_ptr<vector<LorentzVector> >	els_vertex_p4				(new vector<LorentzVector>	) ;

  //HitPattern information
  //
  auto_ptr<vector<LorentzVector> >	els_inner_position			(new vector<LorentzVector>	) ;
  auto_ptr<vector<LorentzVector> >	els_outer_position			(new vector<LorentzVector>	) ;
  auto_ptr<vector<int> >		els_valid_pixelhits			(new vector<int>	        ) ; 
  auto_ptr<vector<int> >		els_lost_pixelhits			(new vector<int>        	) ; 
  auto_ptr<vector<int> >		els_layer1_sizerphi			(new vector<int>	        ) ; 
  auto_ptr<vector<int> >		els_layer1_sizerz			(new vector<int>	        ) ; 
  auto_ptr<vector<float> >		els_layer1_charge			(new vector<float>		) ;
  auto_ptr<vector<int> >		els_layer1_det				(new vector<int>	        ) ;
  auto_ptr<vector<int> >		els_layer1_layer			(new vector<int>		) ;
  auto_ptr<vector<int> >		els_exp_innerlayers			(new vector<int>		) ; 
  auto_ptr<vector<int> >		els_exp_outerlayers			(new vector<int>		) ; 

  auto_ptr<vector<int>     >		els_trkidx				(new vector<int>		) ;
  auto_ptr<vector<float>   >		els_trkshFrac				(new vector<float>		) ;
  auto_ptr<vector<float>   >		els_trkdr				(new vector<float>		) ;

  //conversions
  auto_ptr<vector<vector<float> > >		els_convs_dist			( new vector<vector<float> >		) ;
  auto_ptr<vector<vector<float> > >		els_convs_dcot			( new vector<vector<float> >		) ;
  auto_ptr<vector<vector<float> > >		els_convs_radius                ( new vector<vector<float> >		) ;
  auto_ptr<vector<vector<LorentzVector> > >	els_convs_pos_p4		( new vector<vector<LorentzVector> >	) ;
  auto_ptr<vector<vector<int> > >          	els_convs_tkidx			( new vector<vector<int> >		) ;
  auto_ptr<vector<vector<int> > >		els_convs_gsftkidx		( new vector<vector<int> >		) ;
  auto_ptr<vector<vector<int> > >		els_convs_delMissHits		( new vector<vector<int> >		) ;
  auto_ptr<vector<vector<int> > >         	els_convs_flag			( new vector<vector<int> >		) ;  
      
  auto_ptr< vector <float>         >	els_conv_dist			( new vector<float>		) ;
  auto_ptr< vector <float>         >	els_conv_dcot			( new vector<float>		) ;
  auto_ptr< vector <float>         >	els_conv_radius			( new vector<float>		) ;
  auto_ptr< vector <LorentzVector> >	els_conv_pos_p4			( new vector<LorentzVector>	) ;
  auto_ptr< vector <int>           >	els_conv_tkidx			( new vector<int>		) ;
  auto_ptr< vector <int>           >	els_conv_gsftkidx		( new vector<int>		) ;
  auto_ptr< vector <int>           >	els_conv_delMissHits		( new vector<int>		) ;
  auto_ptr< vector <int>           >	els_conv_flag			( new vector<int>		) ;

  auto_ptr< vector <float>         >	els_conv_old_dist		( new vector<float>		) ;
  auto_ptr< vector <float>         >	els_conv_old_dcot		( new vector<float>		) ;
  auto_ptr< vector <float>         >	els_conv_old_radius		( new vector<float>		) ;
  auto_ptr< vector <int>           >	els_conv_old_tkidx		( new vector<int>		) ;
  auto_ptr< vector <int>           >	els_conv_old_gsftkidx		( new vector<int>		) ;
  auto_ptr< vector <int>           >	els_conv_old_delMissHits	( new vector<int>		) ;
  auto_ptr< vector <int>           >	els_conv_old_flag		( new vector<int>		) ;

  

   
  Handle<reco::TrackCollection> tracks_h;
  iEvent.getByLabel(trksInputTag_, tracks_h);

  //get GSF Tracks
  Handle<reco::GsfTrackCollection> gsftracks_h;
  iEvent.getByLabel(gsftracksInputTag_, gsftracks_h);

  Handle<float> evt_bField_h;
  iEvent.getByLabel("eventMaker", "evtbField", evt_bField_h);
  float evt_bField = *evt_bField_h.product();
  
  // Get products from the reco
  //

  // Get the electrons
  //
  Handle<View<reco::GsfElectron> > els_h;
  iEvent.getByLabel(electronsInputTag_, els_h);
  View<reco::GsfElectron> gsfElColl = *(els_h.product());

  // Get tools to get cluster shape information
  //
  if (clusterTools_) delete clusterTools_;
  clusterTools_ = new EcalClusterLazyTools(iEvent, iSetup, 
					   edm::InputTag("reducedEcalRecHitsEB"), 
					   edm::InputTag("reducedEcalRecHitsEE"));

  // Get beamspot
  //
  edm::InputTag beamSpot_tag(beamSpotInputTag_.label(),"evtbsp4");
  edm::Handle<LorentzVector> beamSpotH;
  iEvent.getByLabel(beamSpot_tag, beamSpotH);
  const Point beamSpot = beamSpotH.isValid() ?
    Point(beamSpotH->x(), beamSpotH->y(), beamSpotH->z()) : Point(0,0,0);

  // Get the value maps for the Egamma electron ID decisions
  /*
  const edm::ValueMap<float>&  eidRobustLooseMap          = getValueMap<float>(iEvent, eidRobustLooseTag_);
  const edm::ValueMap<float>&  eidRobustTightMap          = getValueMap<float>(iEvent, eidRobustTightTag_);
  const edm::ValueMap<float>&  eidRobustHighEnergyMap     = getValueMap<float>(iEvent, eidRobustHighEnergyTag_);
  const edm::ValueMap<float>&  eidLooseMap                = getValueMap<float>(iEvent, eidLooseTag_);
  const edm::ValueMap<float>&  eidTightMap                = getValueMap<float>(iEvent, eidTightTag_);
  */
  //get cms2scsseeddetid 
  edm::InputTag cms2scsseeddetid_tag(cms2scsseeddetidInputTag_.label(),"scsdetIdSeed");
  edm::Handle<std::vector<int> > cms2scsseeddetid_h;
  iEvent.getByLabel(cms2scsseeddetid_tag, cms2scsseeddetid_h); 
  const vector<int> *cms2scsseeddetid = cms2scsseeddetid_h.product();


  //fill number of eqlectrons variable
  //
  *evt_nels = els_h->size();

  //loop over electron collection
  //
  size_t elsIndex = 0;
  View<reco::GsfElectron>::const_iterator el;
  for(el = els_h->begin(); el != els_h->end(); el++, elsIndex++) {

    // Get electron and track objects
    const reco::Track *el_track = (const reco::Track*)(el->gsfTrack().get());
    const edm::RefToBase<reco::GsfElectron> gsfElRef = els_h->refAt(elsIndex);    

    // Get cluster info that is not stored in the object
    //
    float eMax = -9999.;
    float e3x3 = -9999.;    
    if(el->superCluster()->seed().isAvailable() ) { 
      const reco::BasicCluster& clRef= *(el->superCluster()->seed());
      eMax = clusterTools_->eMax(clRef);
      e3x3 = clusterTools_->e3x3(clRef);
    

      // get the covariances computed in 5x5 around the seed
      const std::vector<float>& covs = clusterTools_->covariances(clRef);
      // get the local covariances computed in a 5x5 around the seed
      const std::vector<float>& lcovs = clusterTools_->localCovariances(clRef);
      // get the local covariances computed using all crystals in the SC
      const std::vector<float> localCovariancesSC = clusterTools_->scLocalCovariances(*(el->superCluster()));

      els_eSeed          ->push_back( el->superCluster()->seed()->energy()     );		
      els_sigmaPhiPhi    ->push_back( std::isfinite(covs[2])               ? covs[2] > 0                ? sqrt(covs[2])  : -1 * sqrt(-1 * covs[2])                              : -9999. );
      els_sigmaIPhiIPhi  ->push_back( std::isfinite(lcovs[2])              ? lcovs[2] > 0               ? sqrt(lcovs[2]) : -1 * sqrt(-1 * lcovs[2])                             : -9999. );
      els_sigmaIEtaIEtaSC->push_back( std::isfinite(localCovariancesSC[0]) ? localCovariancesSC[0] > 0  ? sqrt(localCovariancesSC[0])   : -1 * sqrt(-1 * localCovariancesSC[0]) : -9999. );
      els_sigmaIPhiIPhiSC->push_back( std::isfinite(localCovariancesSC[2]) ? localCovariancesSC[2] > 0  ? sqrt(localCovariancesSC[2])   : -1 * sqrt(-1 * localCovariancesSC[2]) : -9999. );
    } else {
      els_eSeed			->push_back(	-9999.	);
      els_sigmaPhiPhi		->push_back(	-9999.	);
      els_sigmaIPhiIPhi		->push_back(	-9999.	);
      els_sigmaIEtaIEtaSC	->push_back(	-9999.	);
      els_sigmaIPhiIPhiSC	->push_back(	-9999.	);
    }
    

    
    // Fill cluster info
    //
    els_etaSC		   ->push_back( el->superCluster()->eta()		);
    els_phiSC		   ->push_back( el->superCluster()->phi()	       	);
    els_eSC            ->push_back( el->superCluster()->energy()             );
    els_eSCRaw         ->push_back( el->superCluster()->rawEnergy()          );
    els_eSCPresh       ->push_back( el->superCluster()->preshowerEnergy()    );
    els_nSeed          ->push_back( el->basicClustersSize() - 1              );      
    els_e1x5		   ->push_back( el->e1x5()				);
    els_e3x3           ->push_back( e3x3                                     );
    els_e5x5           ->push_back( el->e5x5()                               );
    els_e2x5Max        ->push_back( el->e2x5Max()                            );
    els_eMax           ->push_back( eMax                                     );
    
    els_sigmaEtaEta    ->push_back( el->scSigmaEtaEta()                      );
    els_sigmaIEtaIEta  ->push_back( el->scSigmaIEtaIEta()                    );


    //get cms2scsseeddetid--sorry for junk from photons...
    bool foundseed = false;
    for( unsigned int i=0; i<cms2scsseeddetid->size(); i++ ) {      
      if( cms2scsseeddetid->at(i) == -9999)
	continue;      
      if(!(el->superCluster()->seed().isAvailable()) )
	continue;
      if(  uint32_t(cms2scsseeddetid->at(i)) == el->superCluster()->seed()->seed() ) {
	foundseed = true;
	els_scindex->push_back( i );
	break;
      }
    }
    //cout << endl;
    if( !foundseed ) {
      //this is understood: the photon can have energy significantly higher than SC for whatever reason.
      //cout << "No seed found. seed id: " << int(photon->superCluster()->seed()->seed())
      //	 << "  photon et: " << photon->et()
      //	 << "  sc et: " << photon->superCluster()->energy()/cosh(photon->superCluster()->eta())
	  //	 << endl;
      els_scindex->push_back( -1 );		
    }

    // set the mask that describes the egamma fiduciality flags
    // the enum is in interface/EgammaFiduciality.h
    int fiducialityMask = 0;
    if (el->isEB())     	fiducialityMask |= 1 << ISEB;
    if (el->isEBEEGap())	fiducialityMask |= 1 << ISEBEEGAP;
    if (el->isEE())     	fiducialityMask |= 1 << ISEE;
    if (el->isEEGap())  	fiducialityMask |= 1 << ISEEGAP;
    if (el->isEBEtaGap()) 	fiducialityMask |= 1 << ISEBETAGAP;
    if (el->isEBPhiGap())	fiducialityMask |= 1 << ISEBPHIGAP;
    if (el->isEEDeeGap())	fiducialityMask |= 1 << ISEEDEEGAP;
    if (el->isEERingGap())	fiducialityMask |= 1 << ISEERINGGAP;
    if (el->isGap())	fiducialityMask |= 1 << ISGAP;
    els_fiduciality->push_back( fiducialityMask );
			    
    // what corrections have been applied
    // and how is the electron seeding driven
    int electronTypeMask = 0;
    if (el->isEcalEnergyCorrected()) 	electronTypeMask |= 1 << ISECALENERGYCORRECTED;
    // Depricated in CMSSW_4_2x ( http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DataFormats/EgammaCandidates/interface/GsfElectron.h?revision=1.50&view=markup )
    //if (el->isMomentumCorrected())		electronTypeMask |= 1 << ISMOMENTUMCORRECTED; 
    if (el->trackerDrivenSeed())		electronTypeMask |= 1 << ISTRACKERDRIVEN;
    if (el->ecalDrivenSeed())			electronTypeMask |= 1 << ISECALDRIVEN; 
    els_type->push_back( electronTypeMask);

    // energy corrections and uncertainties
    els_ecalEnergy		->push_back(	el->ecalEnergy()	        	);
    els_ecalEnergyError		->push_back(	el->ecalEnergyError()		        );
    els_trackMomentumError	->push_back(	el->trackMomentumError()	        );
    //els_electronMomentumError	->push_back(	el->electronMomentumError()	        );

    // Fill predifined electron ID decisions
    //
    // this is the old pTDR classification
    els_class                 ->push_back( el->classification()	            	); 
    //
    // this is the sani classification
    els_category              ->push_back( classify(gsfElRef)		        );
    //
    // Now get the decisions from the "official" Egamma sequences

    // Track parameters
    //
    float pt = el_track->pt();
    float p = el_track->p();
    float q = el_track->charge();
    float pz = el_track->pz();
    float trkpterr = (el_track->charge()!=0) ? sqrt(pt*pt*p*p/pow(q, 2)*(el_track->covariance(0,0))
						    +2*pt*p/q*pz*(el_track->covariance(0,1))
						    + pz*pz*(el_track->covariance(1,1) ) ) : -9999.;
            
    els_chi2                  ->push_back( el_track->chi2()                          );
    els_ndof                  ->push_back( el_track->ndof()                          );
    els_d0Err                 ->push_back( el_track->d0Error()                       );
    els_z0Err                 ->push_back( el_track->dzError()                       );
    els_ptErr                 ->push_back( trkpterr                                  );
    els_etaErr                ->push_back( el_track->etaError()                      );
    els_phiErr                ->push_back( el_track->phiError()                      );  		
    els_gsftrkidx             ->push_back( static_cast<int>((el->gsfTrack()).key())  );


    els_validHits             ->push_back( el_track->numberOfValidHits()             );
    els_lostHits              ->push_back( el_track->numberOfLostHits()              );
    
    els_charge                ->push_back( el->charge()                              );
    els_trk_charge            ->push_back( el_track->charge()                        );
    els_sccharge              ->push_back( el->scPixCharge()                         );
    els_d0                    ->push_back( el_track->d0()                            );
    els_z0                    ->push_back( el_track->dz()                            );		
    els_d0corr                ->push_back( -1*(el_track->dxy(beamSpot))              );
    els_z0corr                ->push_back( el_track->dz(beamSpot)                    );
    
    // Lorentz Vectors	
    //
    LorentzVector trk_p4( el_track->px(), el_track->py(), 
			  el_track->pz(), el_track->p() );
    double          mass = 0.000510998918;
    LorentzVector   p4In; 
    math::XYZVectorF p3In = el->trackMomentumAtVtx();
    p4In.SetXYZT(   p3In.x(), p3In.y(), p3In.z(), sqrt(mass*mass+p3In.R()*p3In.R()));
    LorentzVector   p4Out; 
    math::XYZVectorF p3Out = el->trackMomentumOut();
    p4Out.SetXYZT(  p3Out.x(), p3Out.y(), p3Out.z(), sqrt(mass*mass+p3Out.R()*p3Out.R()));

    els_p4                    ->push_back( LorentzVector( el->p4() )                 );
    els_trk_p4                ->push_back( trk_p4                                    );
    els_p4In                  ->push_back( p4In                                      );
    els_p4Out                 ->push_back( p4Out                                     );

    // Isolation related
    //
    els_ecalIso               ->push_back(el->dr03EcalRecHitSumEt()                  );
    els_hcalIso               ->push_back(el->dr03HcalTowerSumEt()                   );
    els_hcalDepth1TowerSumEt  ->push_back(el->dr03HcalDepth1TowerSumEt()	     );
    els_hcalDepth2TowerSumEt  ->push_back(el->dr03HcalDepth2TowerSumEt()             );
    els_tkIso                 ->push_back(el->dr03TkSumPt()                          );

    els_ecalIso04             ->push_back(el->dr04EcalRecHitSumEt()                  );
    els_hcalIso04             ->push_back(el->dr04HcalTowerSumEt()                   );
    els_hcalDepth1TowerSumEt04->push_back(el->dr04HcalDepth1TowerSumEt()             );
    els_hcalDepth2TowerSumEt04->push_back(el->dr04HcalDepth2TowerSumEt()             );
    els_tkIso04               ->push_back(el->dr04TkSumPt()                          );

    // Electron ID variables
    //
    double phi_pin = el->caloPosition().phi() - el->deltaPhiSuperClusterTrackAtVtx();
    double phi_pout = -9999.;
    if(el->superCluster()->seed().isAvailable())
      phi_pout = el->superCluster()->seed()->position().phi() - el->deltaPhiSeedClusterTrackAtCalo();
    els_hOverE                ->push_back( el->hadronicOverEm()                      );
    els_hcalDepth1OverEcal    ->push_back( el->hcalDepth1OverEcal()		     );
    els_hcalDepth2OverEcal    ->push_back( el->hcalDepth2OverEcal()                  );

    els_eOverPIn              ->push_back( el->eSuperClusterOverP()                  );
    els_eSeedOverPOut         ->push_back( el->eSeedClusterOverPout()                );
    els_eSeedOverPIn          ->push_back( el->eSeedClusterOverP()		     );
    els_eOverPOut  	      ->push_back( el->eEleClusterOverPout()		     );

    els_fbrem                 ->push_back( el->fbrem()                            );

    els_mva                     ->push_back( el->mva()                                  );
    els_dEtaIn	                ->push_back( el->deltaEtaSuperClusterTrackAtVtx()	    );
    els_dEtaOut                 ->push_back( el->deltaEtaSeedClusterTrackAtCalo()       );
    els_deltaEtaEleClusterTrackAtCalo ->push_back(el->deltaEtaEleClusterTrackAtCalo()   );
    els_dPhiIn                ->push_back( el->deltaPhiSuperClusterTrackAtVtx()         );		
    els_dPhiInPhiOut          ->push_back( phi_pin - phi_pout                           );
    els_dPhiOut               ->push_back( el->deltaPhiSeedClusterTrackAtCalo()         );
    els_deltaPhiEleClusterTrackAtCalo ->push_back( el->deltaPhiEleClusterTrackAtCalo()  );

    // Vertex
    //
    els_vertex_p4             ->push_back( LorentzVector(el->vx(), el->vy(), el->vz(), 0.) );

    //Hit Pattern 
    if(el_track->extra().isAvailable()) {
      els_inner_position ->push_back(LorentzVector(el_track->innerPosition().x(), el_track->innerPosition().y() , el_track->innerPosition().z(), 0 ));
      els_outer_position ->push_back(LorentzVector(el_track->outerPosition().x(), el_track->outerPosition().y() , el_track->outerPosition().z(), 0 ));
    } else {
      els_inner_position->push_back(LorentzVector(-9999., -9999., -9999., -9999.));
      els_outer_position->push_back(LorentzVector(-9999., -9999., -9999., -9999.));
    }
    
    const reco::HitPattern& pattern = el_track->hitPattern();
    const reco::HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 
    const reco::HitPattern& p_outer = el_track->trackerExpectedHitsOuter();
    els_exp_innerlayers    -> push_back(p_inner.numberOfHits());
    els_exp_outerlayers    -> push_back(p_outer.numberOfHits());
    els_valid_pixelhits ->push_back(pattern.numberOfValidPixelHits());
    els_lost_pixelhits ->push_back(pattern.numberOfLostPixelHits());

    if(el_track->extra().isAvailable()) {
      bool valid_hit      = false;
      uint32_t hit_pattern; 
      int i_layer       = 1;
      int side = -1;
      bool pixel_hit   = false;
      bool strip_hit   = false;

      int pixel_size;
      int pixel_sizeX;
      int pixel_sizeY;
      float pixel_charge;
      int det;
      int layer;

      typedef edm::Ref<edmNew::DetSetVector<SiStripCluster>,SiStripCluster > ClusterRef;
      typedef edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster > pixel_ClusterRef;


      for(trackingRecHit_iterator ihit = el_track->recHitsBegin(); 
	  ihit != el_track->recHitsEnd(); ++ihit){
	if(i_layer > 1) break;
	int k = ihit-el_track->recHitsBegin();
	hit_pattern = pattern.getHitPattern(k);
	valid_hit = pattern.validHitFilter(hit_pattern);
	pixel_hit = pattern.pixelHitFilter(hit_pattern);
	strip_hit = pattern.stripHitFilter(hit_pattern);
	side      = (int)pattern.getSide(hit_pattern);
	det       = (int)pattern.getSubStructure(hit_pattern);
	layer     = (int)pattern.getLayer(hit_pattern);

	if(!valid_hit) continue;
	if(pixel_hit){

	  const SiPixelRecHit *pixel_hit_cast = dynamic_cast<const SiPixelRecHit*>(&(**ihit));
	  assert(pixel_hit_cast != 0);
	  pixel_ClusterRef const& pixel_cluster = pixel_hit_cast->cluster();

	  pixel_size   = (int)pixel_cluster->size(); 
	  pixel_sizeX  = (int)pixel_cluster->sizeX(); 
	  pixel_sizeY  = (int)pixel_cluster->sizeY(); 
	  pixel_charge = (float)pixel_cluster->charge();

	  if(i_layer == 1){
	    els_layer1_sizerphi ->push_back(pixel_sizeX);
	    els_layer1_sizerz   ->push_back(pixel_sizeY);
	    els_layer1_charge   ->push_back(pixel_charge);
	    els_layer1_det      ->push_back(det);
	    els_layer1_layer    ->push_back(layer);
	    i_layer++;

	  }
	}

	else if (strip_hit){
	  const SiStripRecHit1D *strip_hit_cast = dynamic_cast<const SiStripRecHit1D*>(&(**ihit));
	  const SiStripRecHit2D *strip2d_hit_cast = dynamic_cast<const SiStripRecHit2D*>(&(**ihit));
	  ClusterRef cluster;
	  if(strip_hit_cast == NULL)
	    cluster = strip2d_hit_cast->cluster();
	  else 
	    cluster = strip_hit_cast->cluster();

	  int cluster_size   = (int)cluster->amplitudes().size();
	  int cluster_charge = 0;
	  double   cluster_weight_size = 0.0;
	  int max_strip_i = std::max_element(cluster->amplitudes().begin(),cluster->amplitudes().end())-cluster->amplitudes().begin();

	  for(int istrip = 0; istrip < cluster_size; istrip++){
	    cluster_charge += (int)cluster->amplitudes().at(istrip);
	    cluster_weight_size += (istrip-max_strip_i)*(istrip-max_strip_i)*(cluster->amplitudes().at(istrip));
	  }
	  cluster_weight_size = sqrt(cluster_weight_size/cluster_charge);

	  if(i_layer == 1){
	    if(side==0) 
	      {
		els_layer1_sizerphi ->push_back(cluster_size);
		els_layer1_sizerz   ->push_back(0);
	      }

	    else
	      {
		els_layer1_sizerphi ->push_back(0);
		els_layer1_sizerz   ->push_back(cluster_size);
	      } 

	    els_layer1_charge   ->push_back(cluster_charge);
	    els_layer1_det      ->push_back(det);
	    els_layer1_layer    ->push_back(layer);
	    i_layer++;
	  }
	}

      }


    } else {
      els_layer1_sizerphi ->push_back(-9999);
      els_layer1_sizerz   ->push_back(-9999);
      els_layer1_charge   ->push_back(-9999);
      els_layer1_det      ->push_back(-9999);
      els_layer1_layer    ->push_back(-9999);
    }
    
    // *****************************************************

    TrackRef ctfTkRef    = el->closestCtfTrackRef();
    GsfTrackRef gsfTkRef = el->gsfTrack();

    double dR = -9999.;
    if(ctfTkRef.isNonnull() ) {
      els_trkidx       ->push_back(static_cast<int>(ctfTkRef.key())            );
      els_trkshFrac    ->push_back(static_cast<float>(el->shFracInnerHits())     );
      dR = deltaR(gsfTkRef->eta(), gsfTkRef->phi(),
		  ctfTkRef->eta(), ctfTkRef->phi()                             );
    } else {
      els_trkidx       ->push_back(-9999                                        );
      els_trkshFrac    ->push_back(-9999.                                       );
    }

    els_trkdr          ->push_back(dR                                          );




    ConversionFinder convFinder;
    //vector of conversion infos - all the candidate conversion tracks
    vector<ConversionInfo> v_convInfos = convFinder.getConversionInfos(*(el->core()), tracks_h, gsftracks_h, evt_bField);
    vector<float> v_dist, v_dcot, v_rad;
    vector<int>   v_tkidx, v_gsftkidx, v_delmisshits, v_flag;
    vector<LorentzVector> v_pos_p4;

    for(unsigned int i_conv = 0; i_conv < v_convInfos.size(); i_conv++) {
      
      v_dist.		push_back(std::isfinite(v_convInfos.at(i_conv).dist()) ? v_convInfos.at(i_conv).dist() : -9999.);
      v_dcot.		push_back(v_convInfos.at(i_conv).dcot());
      v_rad.		push_back(v_convInfos.at(i_conv).radiusOfConversion());
      math::XYZPoint convPoint		= v_convInfos.at(i_conv).pointOfConversion();
      float convPointX			= std::isfinite(convPoint.x()) ? convPoint.x() : -9999.;
      float convPointY			= std::isfinite(convPoint.y()) ? convPoint.y() : -9999.;
      float convPointZ			= std::isfinite(convPoint.z()) ? convPoint.z() : -9999.;
      v_pos_p4.		push_back(LorentzVector(convPointX, convPointY, convPointZ, 0));
      if(v_convInfos.at(i_conv).conversionPartnerCtfTk().isNonnull())
	v_tkidx.	push_back(v_convInfos.at(i_conv).conversionPartnerCtfTk().key());
      else 
	v_tkidx.	push_back(-9999);
      if(v_convInfos.at(i_conv).conversionPartnerGsfTk().isNonnull())
	v_gsftkidx.	push_back(v_convInfos.at(i_conv).conversionPartnerGsfTk().key());
      else 
	v_gsftkidx.	push_back(-9999);
      v_delmisshits.	push_back(v_convInfos.at(i_conv).deltaMissingHits());
      v_flag.		push_back(v_convInfos.at(i_conv).flag());

    }

    els_convs_dist		->push_back(v_dist		);
    els_convs_dcot		->push_back(v_dcot		);
    els_convs_radius		->push_back(v_rad		);
    els_convs_pos_p4		->push_back(v_pos_p4		);
    els_convs_tkidx		->push_back(v_tkidx		);
    els_convs_gsftkidx		->push_back(v_gsftkidx		);
    els_convs_delMissHits	->push_back(v_delmisshits	);
    els_convs_flag		->push_back(v_flag		);



    ConversionInfo convInfo	= convFinder.getConversionInfo(*el, tracks_h, gsftracks_h, evt_bField);
    els_conv_dist		->push_back(std::isfinite(convInfo.dist()) ? convInfo.dist() : -9999.);
    els_conv_dcot		->push_back(convInfo.dcot());
    els_conv_radius		->push_back(convInfo.radiusOfConversion());
    math::XYZPoint convPoint	= convInfo.pointOfConversion();
    float convPointX		= std::isfinite(convPoint.x()) ? convPoint.x() : -9999.;
    float convPointY		= std::isfinite(convPoint.y()) ? convPoint.y() : -9999.;
    float convPointZ		= std::isfinite(convPoint.z()) ? convPoint.z() : -9999.;
    els_conv_pos_p4		->push_back(LorentzVector(convPointX, convPointY, convPointZ, 0));
    if(convInfo.conversionPartnerCtfTk().isNonnull())
      els_conv_tkidx		->push_back(convInfo.conversionPartnerCtfTk().key());
    else 
      els_conv_tkidx		->push_back(-9999);
    if(convInfo.conversionPartnerGsfTk().isNonnull())
      els_conv_gsftkidx		->push_back(convInfo.conversionPartnerGsfTk().key());
    else 
      els_conv_gsftkidx		->push_back(-9999);
    els_conv_delMissHits	->push_back(convInfo.deltaMissingHits());
    els_conv_flag		->push_back(convInfo.flag());

    // Old Conversion Rejection
    int delMissHits          =  -9999;
    int flag                 =  isfinite(el->convFlags()) ? el->convFlags() : -9999;
    els_conv_old_flag        -> push_back( flag                     );
    els_conv_old_dist        -> push_back( isfinite(el->convDist())   ? el->convDist()   : -9999. );
    els_conv_old_dcot        -> push_back( isfinite(el->convDcot())   ? el->convDcot()   : -9999. );
    els_conv_old_radius      -> push_back( isfinite(el->convRadius()) ? el->convRadius() : -9999. );
    if( flag == 0) { //partner found in the CTF track collection using the electron's CTF track
      els_conv_old_tkidx -> push_back( el->convPartner().key() );
      els_conv_old_gsftkidx-> push_back( -9999 );
      delMissHits = el->convPartner()->trackerExpectedHitsInner().numberOfHits() - el->closestCtfTrackRef()->trackerExpectedHitsInner().numberOfHits();
    }
    else if( flag == 1) {//partner found in the CTF track collection using the electron's GSF track
      els_conv_old_tkidx -> push_back( el->convPartner().key() );  
      els_conv_old_gsftkidx-> push_back( -9999 );
      delMissHits = el->convPartner()->trackerExpectedHitsInner().numberOfHits() - el_track->trackerExpectedHitsInner().numberOfHits();
    }
    else if(flag == 2) {//partner found in the GSF track collection using the electron's CTF track
      els_conv_old_tkidx -> push_back( -9999 );
      els_conv_old_gsftkidx ->push_back(el->convPartner().key() );
      delMissHits = el->convPartner()->trackerExpectedHitsInner().numberOfHits() - el->closestCtfTrackRef()->trackerExpectedHitsInner().numberOfHits();
    } else if(flag == 3) {//partner found in the GSF track collection using the electron's GSF track
      els_conv_old_tkidx -> push_back( -9999 );
      els_conv_old_gsftkidx ->push_back(el->convPartner().key() );
      delMissHits = el->convPartner()->trackerExpectedHitsInner().numberOfHits() - el_track->trackerExpectedHitsInner().numberOfHits();
    } else if( flag != -9999 ){
      cout << "ERROR Unexpected flag: " << flag << endl;
      exit(1);
    }
    els_conv_old_delMissHits -> push_back( delMissHits              );


  }//electron loop
  


  // Put the results into the event
  //
  iEvent.put(evt_nels           	    	,"evtnels"				);

  // Predefined ID descisions 
  //
  iEvent.put(els_class                   	,"elsclass"				);
  iEvent.put(els_category                	,"elscategory"				);
  // Track parameters
  //
  iEvent.put(els_d0                     	,"elsd0"				);
  iEvent.put(els_z0                      	,"elsz0"				);
  iEvent.put(els_d0corr				,"elsd0corr"				);
  iEvent.put(els_z0corr				,"elsz0corr"				);
  iEvent.put(els_chi2                    	,"elschi2"				);
  iEvent.put(els_ndof                   	,"elsndof"				);
  iEvent.put(els_d0Err				,"elsd0Err"				);
  iEvent.put(els_z0Err                   	,"elsz0Err"				);
  iEvent.put(els_ptErr                   	,"elsptErr"				);
  iEvent.put(els_etaErr                  	,"elsetaErr"				);
  iEvent.put(els_phiErr                  	,"elsphiErr"				);
  iEvent.put(els_gsftrkidx                      ,"elsgsftrkidx"                         );
  
  iEvent.put(els_validHits               	,"elsvalidHits"				);
  iEvent.put(els_lostHits                	,"elslostHits"				);
  iEvent.put(els_charge                  	,"elscharge"				);
  iEvent.put(els_trk_charge                     ,"elstrkcharge"				);
  iEvent.put(els_sccharge                       ,"elssccharge"				);

  // Supercluster parameters
  //
  iEvent.put(els_nSeed				,"elsnSeed"				);
  iEvent.put(els_etaSC				,"elsetaSC"				);
  iEvent.put(els_phiSC				,"elsphiSC"				);
  iEvent.put(els_eSC				,"elseSC"				);
  iEvent.put(els_eSCRaw				,"elseSCRaw"				);
  iEvent.put(els_eSCPresh			,"elseSCPresh"				);
  iEvent.put(els_e1x5				,"else1x5"				);
  iEvent.put(els_e3x3				,"else3x3"				);
  iEvent.put(els_e5x5				,"else5x5"				);
  iEvent.put(els_e2x5Max			,"else2x5Max"				);
  iEvent.put(els_eMax				,"elseMax"				);
  iEvent.put(els_eSeed				,"elseSeed"				);
  iEvent.put(els_fiduciality			,"elsfiduciality"			);
  iEvent.put(els_type				,"elstype"				);
  iEvent.put(els_scindex			,"elsscindex"				);

  // Corrections and uncertainties
  //
  iEvent.put(els_ecalEnergy			,"elsecalEnergy"			);
  iEvent.put(els_ecalEnergyError		,"elsecalEnergyError"			);
  iEvent.put(els_trackMomentumError		,"elstrackMomentumError"		);
  iEvent.put(els_electronMomentumError		,"elselectronMomentumError"		);

  // Electron ID
  //
  iEvent.put(els_sigmaPhiPhi			,"elssigmaPhiPhi"			);
  iEvent.put(els_sigmaIPhiIPhi			,"elssigmaIPhiIPhi"			);
  iEvent.put(els_sigmaEtaEta			,"elssigmaEtaEta"			);
  iEvent.put(els_sigmaIEtaIEta			,"elssigmaIEtaIEta"			);
  iEvent.put(els_sigmaIPhiIPhiSC		,"elssigmaIPhiIPhiSC"			);
  iEvent.put(els_sigmaIEtaIEtaSC		,"elssigmaIEtaIEtaSC"			);
  iEvent.put(els_dPhiInPhiOut			,"elsdPhiInPhiOut"			);
  iEvent.put(els_hOverE				,"elshOverE"				);
  iEvent.put(els_hcalDepth1OverEcal		,"elshcalDepth1OverEcal"		);
  iEvent.put(els_hcalDepth2OverEcal		,"elshcalDepth2OverEcal"		);

  iEvent.put(els_eOverPIn			,"elseOverPIn"				);
  iEvent.put(els_eSeedOverPOut			,"elseSeedOverPOut"			);
  iEvent.put(els_eSeedOverPIn			,"elseSeedOverPIn"			);
  iEvent.put(els_eOverPOut			,"elseOverPOut"				);
  iEvent.put(els_fbrem				,"elsfbrem"				);
  iEvent.put(els_mva				,"elsmva"				);
  iEvent.put(els_dEtaIn				,"elsdEtaIn"				);
  iEvent.put(els_dEtaOut			,"elsdEtaOut"				);
  iEvent.put(els_deltaEtaEleClusterTrackAtCalo	, "elsdeltaEtaEleClusterTrackAtCalo"	);
  iEvent.put(els_dPhiIn				,"elsdPhiIn"				);
  iEvent.put(els_dPhiOut			,"elsdPhiOut"				);
  iEvent.put(els_deltaPhiEleClusterTrackAtCalo	, "elsdeltaPhiEleClusterTrackAtCalo"	);

  // Lorentz vectors
  //
  iEvent.put(els_p4				,"elsp4"				);
  iEvent.put(els_trk_p4				,"elstrkp4"				);
  iEvent.put(els_p4In				,"elsp4In"				);
  iEvent.put(els_p4Out				,"elsp4Out"				);

  // Vertex
  //
  iEvent.put(els_vertex_p4			,"elsvertexp4"				);

  // Isolation
  //
  iEvent.put(els_tkIso				,"elstkIso"				);
  iEvent.put(els_ecalIso			,"elsecalIso"				);
  iEvent.put(els_hcalIso			,"elshcalIso"				);
  iEvent.put(els_hcalDepth1TowerSumEt		,"elshcalDepth1TowerSumEt"		);
  iEvent.put(els_hcalDepth2TowerSumEt		,"elshcalDepth2TowerSumEt"		);

  iEvent.put(els_tkIso04			,"elstkIso04"				);
  iEvent.put(els_ecalIso04			,"elsecalIso04"				);
  iEvent.put(els_hcalIso04			,"elshcalIso04"				);
  iEvent.put(els_hcalDepth1TowerSumEt04		,"elshcalDepth1TowerSumEt04"		);
  iEvent.put(els_hcalDepth2TowerSumEt04		,"elshcalDepth2TowerSumEt04"		);

  //Hit Pattern Information

  iEvent.put(els_inner_position			,"elsinnerposition"			);
  iEvent.put(els_outer_position			,"elsouterposition"			);
  iEvent.put(els_valid_pixelhits		,"elsvalidpixelhits"			);
  iEvent.put(els_lost_pixelhits			,"elslostpixelhits"			);
  iEvent.put(els_layer1_layer			,"elslayer1layer"			);
  iEvent.put(els_layer1_sizerphi		,"elslayer1sizerphi"			);
  iEvent.put(els_layer1_sizerz			,"elslayer1sizerz"			);
  iEvent.put(els_layer1_charge			,"elslayer1charge"			);
  iEvent.put(els_layer1_det			,"elslayer1det"				);
  iEvent.put(els_exp_innerlayers		,"elsexpinnerlayers"			);
  iEvent.put(els_exp_outerlayers		,"elsexpouterlayers"			);

  //CTF track info
  //
  iEvent.put(els_trkidx                         ,"elstrkidx"				);
  iEvent.put(els_trkdr                          ,"elstrkdr"				);
  iEvent.put(els_trkshFrac                      ,"elstrkshFrac"				);


  //conversion
  iEvent.put(els_convs_dist			, "elsconvsdist"			);
  iEvent.put(els_convs_dcot			, "elsconvsdcot"			);
  iEvent.put(els_convs_radius			, "elsconvsradius"			);
  iEvent.put(els_convs_pos_p4			, "elsconvsposp4"			);
  iEvent.put(els_convs_tkidx			, "elsconvstkidx"			);
  iEvent.put(els_convs_gsftkidx			, "elsconvsgsftkidx"			);
  iEvent.put(els_convs_delMissHits		, "elsconvsdelMissHits"			);
  iEvent.put(els_convs_flag			, "elsconvsflag"			);


  iEvent.put(els_conv_dist			, "elsconvdist"				);
  iEvent.put(els_conv_dcot			, "elsconvdcot"				);
  iEvent.put(els_conv_radius			, "elsconvradius"			);
  iEvent.put(els_conv_pos_p4			, "elsconvposp4"			);
  iEvent.put(els_conv_tkidx			, "elsconvtkidx"			);
  iEvent.put(els_conv_gsftkidx			, "elsconvgsftkidx"			);
  iEvent.put(els_conv_delMissHits		, "elsconvdelMissHits"			);
  iEvent.put(els_conv_flag			, "elsconvflag"				);
  
  iEvent.put(els_conv_old_dist			, "elsconvolddist"			);
  iEvent.put(els_conv_old_dcot			, "elsconvolddcot"			);
  iEvent.put(els_conv_old_radius		, "elsconvoldradius"			);
  iEvent.put(els_conv_old_tkidx			, "elsconvoldtkidx"			);
  iEvent.put(els_conv_old_gsftkidx		, "elsconvoldgsftkidx"			);
  iEvent.put(els_conv_old_delMissHits		, "elsconvolddelMissHits"		);
  iEvent.put(els_conv_old_flag			, "elsconvoldflag"			);

}

//----------------------------------------------------------------------------
// Electron Id classification function (a flag for the Sani type class)
//----------------------------------------------------------------------------
int ElectronMaker::classify(const edm::RefToBase<reco::GsfElectron> &electron) {

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
template<typename T> const edm::ValueMap<T>& ElectronMaker::getValueMap(const edm::Event& iEvent, edm::InputTag& inputTag)
{
  edm::Handle<edm::ValueMap<T> > handle;
  iEvent.getByLabel(inputTag,handle);
  return *(handle.product());
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMaker);

