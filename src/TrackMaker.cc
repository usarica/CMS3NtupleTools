// -*- C++ -*-
//
// Package:    TrackMaker
// Class:      TrackMaker
// 
/**\class TrackMaker TrackMaker.cc CMS2/TrackMaker/src/TrackMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackMaker.cc,v 1.15 2009/08/30 16:29:47 fgolf Exp $
//
//

// system include files
#include <memory>
#include "Math/VectorUtil.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/TrackMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"

// Propagator specific include files
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;
using std::vector;

//
// class decleration
//

//
// constructors and destructor
//
TrackMaker::TrackMaker(const edm::ParameterSet& iConfig)
{
  produces<vector<LorentzVector> >	("trkstrkp4"	  ).setBranchAlias("trks_trk_p4"     );	// track p4						
  produces<vector<LorentzVector> >	("trksvertexp4"	  ).setBranchAlias("trks_vertex_p4"  );	// track p4
  produces<vector<LorentzVector> >      ("trksouterp4"    ).setBranchAlias("trks_outer_p4"   );    // p4 at the outermost point of the tracker
  produces<vector<float> >		("trksd0"	  ).setBranchAlias("trks_d0"         );	// impact parameter at the point of closest approach	
  produces<vector<float> >		("trksd0corr"	  ).setBranchAlias("trks_d0corr"     );	// impact parameter at the point of closest approach corrected for the beamSpot
  produces<vector<float> >		("trksd0corrPhi"  ).setBranchAlias("trks_d0corrPhi"  );	// angle of impact parameter corrected for beamSpot
  produces<vector<float> >		("trksz0"	  ).setBranchAlias("trks_z0"         );	// z position of the point of closest approach		
  produces<vector<float> >		("trksz0corr"	  ).setBranchAlias("trks_z0corr"     );	// z position of the point of closest approach corrected for the the beamSpot
  produces<vector<float> >		("trksvertexphi"  ).setBranchAlias("trks_vertexphi"  );	// phi angle of the point of closest approach		
  produces<vector<float> >		("trkschi2"	  ).setBranchAlias("trks_chi2"       );	// chi2 of the silicon tracker fit			
  produces<vector<float> >		("trksndof"	  ).setBranchAlias("trks_ndof"       );	// number of degrees of freedom of the fit		
  produces<vector<int> >		("trksvalidHits"  ).setBranchAlias("trks_validHits"  );	// number of used hits in the fit			
  produces<vector<int> >		("trkslostHits"	  ).setBranchAlias("trks_lostHits"   );	// number of lost hits in the fit			
  produces<vector<float> >		("trksd0Err"	  ).setBranchAlias("trks_d0Err"      );	// error on the impact parameter			
  produces<vector<float> >		("trksz0Err"	  ).setBranchAlias("trks_z0Err"      );	// error on z position of the point of closest approach	
  produces<vector<float> >		("trksptErr"	  ).setBranchAlias("trks_ptErr"      );	// track Pt error					
  produces<vector<float> >		("trksetaErr"	  ).setBranchAlias("trks_etaErr"     );	// track eta error					
  produces<vector<float> >		("trksphiErr"	  ).setBranchAlias("trks_phiErr"     );	// track phi error					
  produces<vector<int> >		("trkscharge"	  ).setBranchAlias("trks_charge"     );	// charge						
  produces<vector<float> >		("trkstkIso"	  ).setBranchAlias("trks_tkIso"      );	// track isolation like els_tkIso
  produces<vector<int> >                ("trksqualityMask").setBranchAlias("trks_qualityMask"); // mask of quality flags

  tracksInputTag = iConfig.getParameter<edm::InputTag>("tracksInputTag");
  beamSpotTag    = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");

  dRConeMin_   = iConfig.getParameter<double>("trkIsolationdRConeMin"  );
  dRConeMax_   = iConfig.getParameter<double>("trkIsolationdRConeMax"  );
  vtxDiffZMax_ = iConfig.getParameter<double>("trkIsolationVtxDiffZMax");
  tkVtxDMax_   = iConfig.getParameter<double>("trkIsolationTkVtxDMax"  );
  ptMin_       = iConfig.getParameter<double>("trkIsolationPtMin"      );
  nHits_       = iConfig.getParameter<int>   ("trkIsolationNHits"      );

}

void TrackMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  std::auto_ptr<vector<LorentzVector> >	vector_trks_trk_p4	(new vector<LorentzVector>	);
  std::auto_ptr<vector<LorentzVector> >	vector_trks_vertex_p4	(new vector<LorentzVector>	);
  std::auto_ptr<vector<LorentzVector> > vector_trks_outer_p4    (new vector<LorentzVector>      );
  std::auto_ptr<vector<float> >		vector_trks_d0		(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_d0corr      (new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_d0corrPhi   (new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_z0		(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_z0corr	(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_vertexphi	(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_chi2	(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_ndof	(new vector<float>		);      
  std::auto_ptr<vector<int> >		vector_trks_validHits	(new vector<int>		);        
  std::auto_ptr<vector<int> >		vector_trks_lostHits	(new vector<int>		);        
  std::auto_ptr<vector<float> >		vector_trks_d0Err	(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_z0Err	(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_ptErr	(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_etaErr	(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_phiErr	(new vector<float>		);      
  std::auto_ptr<vector<int> >		vector_trks_charge	(new vector<int>		);        
  std::auto_ptr<vector<float> >		vector_trks_outerPhi	(new vector<float>		);      
  std::auto_ptr<vector<float> >		vector_trks_outerEta	(new vector<float>		);      
  std::auto_ptr<vector<float> >         vector_trks_outerPt     (new vector<float>              );
  std::auto_ptr<vector<float> >		vector_trks_tkIso	(new vector<float>		);
  std::auto_ptr<vector<int> >           vector_trks_qualityMask (new vector<int>                );

  // get tracks
  Handle<edm::View<reco::Track> > track_h;
  iEvent.getByLabel(tracksInputTag, track_h);

  if( !track_h.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve track collection";
    edm::LogInfo("OutputInfo") << " TrackMaker cannot continue...!";
    return;
  }

  edm::ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);

  if( !theMagField.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve the magnetic field";
    edm::LogInfo("OutputInfo") << " TrackMaker cannot continue...!";
    return;
  }

  const MagneticField* bf = theMagField.product();

  edm::Handle<LorentzVector> beamSpotH;
  iEvent.getByLabel(beamSpotTag, beamSpotH);

  const Point beamSpot = beamSpotH.isValid() ? Point(beamSpotH->x(), beamSpotH->y(), beamSpotH->z()) : Point(0, 0, 0);

  edm::View<reco::Track>::const_iterator tracks_end = track_h->end();

  for (edm::View<reco::Track>::const_iterator i = track_h->begin(); i != tracks_end; ++i) {

    vector_trks_trk_p4       ->push_back( LorentzVector( i->px(), i->py(), i->pz(), i->p() )       );
    vector_trks_vertex_p4    ->push_back( LorentzVector(i->vx(),i->vy(), i->vz(), 0.)              );
    vector_trks_d0           ->push_back( i->d0()                                                  );
    vector_trks_z0           ->push_back( i->dz()                                                  );
											           
    double corrd0 = beamSpotH.isValid() ? -1 * ( i->dxy(beamSpot) ) : i->d0();		           
    vector_trks_d0corr       ->push_back( corrd0                                                   );
											           
    double corrd0phi = atan2( -1 * corrd0 * sin( i->phi() ), corrd0 * cos( i->phi() ) );           
    vector_trks_d0corrPhi    ->push_back( corrd0phi                                                );
											           
    double corrz0 = beamSpotH.isValid() ? i->dz(beamSpot) : i->dz();			           
    vector_trks_z0corr       ->push_back( corrz0                                                   );
											           
    vector_trks_vertexphi    ->push_back( atan2( i->vy(), i->vx() )                                );
    vector_trks_chi2         ->push_back( i->chi2()                                                );
    vector_trks_ndof         ->push_back( i->ndof()                                                );
    vector_trks_validHits    ->push_back( i->numberOfValidHits()                                   );
    vector_trks_lostHits     ->push_back( i->numberOfLostHits()                                    );
    vector_trks_d0Err        ->push_back( i->d0Error()                                             );
    vector_trks_z0Err        ->push_back( i->dzError()                                             );
    vector_trks_ptErr        ->push_back( i->ptError()                                             );
    vector_trks_etaErr       ->push_back( i->etaError()                                            );
    vector_trks_phiErr       ->push_back( i->phiError()                                            );
    vector_trks_charge       ->push_back( i->charge()                                              );
    vector_trks_tkIso        ->push_back( calculateTrkIsolation(track_h.product(), *i, beamSpot)   );
    vector_trks_qualityMask  ->push_back( i->qualityMask()                                         );
	  
    GlobalPoint  tpVertex   ( i->vx(), i->vy(), i->vz() );
    GlobalVector tpMomentum ( i->px(), i->py(), i->pz() );
    int tpCharge ( i->charge() );

    FreeTrajectoryState fts ( tpVertex, tpMomentum, tpCharge, bf);

    const float zdist = 314.;

    const float radius = 130.;

    const float corner = 1.479;

    Plane::PlanePointer lendcap = Plane::build( Plane::PositionType (0, 0, -zdist), Plane::RotationType () );
    Plane::PlanePointer rendcap = Plane::build( Plane::PositionType (0, 0,  zdist), Plane::RotationType () );

    Cylinder::CylinderPointer barrel = Cylinder::build( Cylinder::PositionType (0, 0, 0), Cylinder::RotationType (), radius);

    AnalyticalPropagator myAP (bf, alongMomentum, 2*M_PI);
	  
    TrajectoryStateOnSurface tsos;
	  
    /*
      Trajectory State is at intersection of cylinder and track, 
      not state at the last hit on the track fit. Shouldn't matter that much.
      Not sure what happens for loopers. Caveat emptor!
    */
	  
    if( i->eta() < -corner ) {
      tsos = myAP.propagate( fts, *lendcap);
    }
    else if( fabs(i->eta()) < corner ) {
      tsos = myAP.propagate( fts, *barrel);
    }
    else if( i->eta() > corner ) {
      tsos = myAP.propagate( fts, *rendcap);
    }

    if(tsos.isValid()) {
      vector_trks_outer_p4->push_back( LorentzVector( tsos.globalMomentum().x(),
						      tsos.globalMomentum().y(),
						      tsos.globalMomentum().z(),
						      tsos.globalMomentum().mag() ) );
    }
    else {
      vector_trks_outer_p4->push_back( LorentzVector( -999., -999., -999., -999.) );

    }
  }

  // store vectors
  iEvent.put(vector_trks_trk_p4       , "trkstrkp4"             );
  iEvent.put(vector_trks_vertex_p4    , "trksvertexp4"          );
  iEvent.put(vector_trks_outer_p4     , "trksouterp4"           );
  iEvent.put(vector_trks_d0           , "trksd0"                );
  iEvent.put(vector_trks_d0corr       , "trksd0corr"            );
  iEvent.put(vector_trks_d0corrPhi    , "trksd0corrPhi"         );
  iEvent.put(vector_trks_z0           , "trksz0"                );
  iEvent.put(vector_trks_z0corr       , "trksz0corr"            );
  iEvent.put(vector_trks_vertexphi    , "trksvertexphi"         );
  iEvent.put(vector_trks_chi2         , "trkschi2"              );
  iEvent.put(vector_trks_ndof         , "trksndof"              );
  iEvent.put(vector_trks_validHits    , "trksvalidHits"         );
  iEvent.put(vector_trks_lostHits     , "trkslostHits"          );
  iEvent.put(vector_trks_d0Err        , "trksd0Err"             );
  iEvent.put(vector_trks_z0Err        , "trksz0Err"             );
  iEvent.put(vector_trks_ptErr        , "trksptErr"             );
  iEvent.put(vector_trks_etaErr       , "trksetaErr"            );
  iEvent.put(vector_trks_phiErr       , "trksphiErr"            );
  iEvent.put(vector_trks_charge       , "trkscharge"            );
  iEvent.put(vector_trks_tkIso        , "trkstkIso"             );
  iEvent.put(vector_trks_qualityMask  , "trksqualityMask"       );
}

// ------------ method called once each job just before starting event loop  ------------
void TrackMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void TrackMaker::endJob() {
}

double TrackMaker::calculateTrkIsolation(const edm::View<reco::Track> *tracks, const reco::Track &input, const Point& beamSpot) {
  //
  // calculate track isolation following electron isolation definition in ElectronMaker.cc
  //
  // sum up all track.pt around track if track fulfills:
  // dR < 0.3
  // dR > 0.01
  // d0 < 0.1
  // dZ < 0.5
  // pT >= 1.0
  // nHits > 7

  double sumPt = 0;

  edm::View<reco::Track>::const_iterator tracks_end = tracks->end();

  for ( edm::View<reco::Track>::const_iterator i = tracks->begin(); i != tracks_end; ++i) {

    double corrz0 = i->dz(beamSpot);
    if ( corrz0 > tkVtxDMax_ ) continue;

    double pt = i->pt();
    if (  pt < ptMin_ ) continue;

    if ( i->numberOfValidHits() <= nHits_ ) continue;

    LorentzVector loopTrack( i->px(), i->py(), i->pz(), i->p() );
    LorentzVector inputTrack( input.px(), input.py(), input.pz(), input.p() );
    double dR = ROOT::Math::VectorUtil::DeltaR(loopTrack, inputTrack);
    if (dR < dRConeMin_ ) continue;
    if (dR > dRConeMax_ ) continue;

    double inputCorrz0 = input.dz(beamSpot);
    double dZ          = fabs( corrz0 - inputCorrz0 );
    if ( dZ >= vtxDiffZMax_) continue;
    sumPt += pt;
  }
  
  return sumPt;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackMaker);
