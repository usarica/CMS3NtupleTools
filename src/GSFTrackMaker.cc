// -*- C++ -*-
//
// Package:    GSFTrackMaker
// Class:      GSFTrackMaker
// 
/**\class GSFTrackMaker GSFTrackMaker.cc CMS2/GSFTrackMaker/src/GSFTrackMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/

// system include files
#include <memory>
#include "Math/VectorUtil.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/GSFTrackMaker.h"

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

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "CondFormats/AlignmentRecord/interface/GlobalPositionRcd.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using std::vector;
using reco::Track;
using reco::TrackBase;

//
// class decleration
//

//
// constructors and destructor
//
GSFTrackMaker::GSFTrackMaker(const edm::ParameterSet& iConfig) {
       
  produces<vector<LorentzVector> >	("gsftrksp4"		).setBranchAlias("gsftrks_p4"     );	// track p4						
  produces<vector<LorentzVector> >	("gsftrksvertexp4"	).setBranchAlias("gsftrks_vertex_p4"  );	// track p4
  produces<vector<LorentzVector> >	  ("gsftrksouterp4"	).setBranchAlias("gsftrks_outer_p4"   );    // p4 at the outermost point of the tracker

  produces<vector<float> >		("gsftrksd0"		).setBranchAlias("gsftrks_d0"         );	// impact parameter at the point of closest approach	
  produces<vector<float> >		("gsftrksd0corr"	).setBranchAlias("gsftrks_d0corr"     );	// impact parameter at the point of closest approach corrected for the beamSpot
  produces<vector<float> >		("gsftrksd0corrPhi"	).setBranchAlias("gsftrks_d0corrPhi"  );	// angle of impact parameter corrected for beamSpot
  produces<vector<float> >		("gsftrksz0"		).setBranchAlias("gsftrks_z0"         );	// z position of the point of closest approach		
  produces<vector<float> >		("gsftrksz0corr"	).setBranchAlias("gsftrks_z0corr"     );	// z position of the point of closest approach corrected for the the beamSpot
  produces<vector<float> >		("gsftrkschi2"		).setBranchAlias("gsftrks_chi2"       );	// chi2 of the silicon tracker fit			
  produces<vector<float> >		("gsftrksndof"		).setBranchAlias("gsftrks_ndof"       );	// number of degrees of freedom of the fit		
  produces<vector<int> >		("gsftrksvalidHits"	).setBranchAlias("gsftrks_validHits"  );	// number of used hits in the fit			
  produces<vector<int> >		("gsftrkslostHits"	).setBranchAlias("gsftrks_lostHits"   );	// number of lost hits in the fit			
  produces<vector<float> >		("gsftrksd0Err"		).setBranchAlias("gsftrks_d0Err"      );	// error on the impact parameter			
  produces<vector<float> >		("gsftrksz0Err"		).setBranchAlias("gsftrks_z0Err"      );	// error on z position of the point of closest approach	
  produces<vector<float> >		("gsftrksptErr"		).setBranchAlias("gsftrks_ptErr"      );	// track Pt error					
  produces<vector<float> >		("gsftrksetaErr"	).setBranchAlias("gsftrks_etaErr"     );	// track eta error					
  produces<vector<float> >		("gsftrksphiErr"	).setBranchAlias("gsftrks_phiErr"     );	// track phi error					
  produces<vector<float> >		("gsftrksd0phiCov"	).setBranchAlias("gsftrks_d0phiCov"   ); // track cov(d0, phi) 
  produces<vector<int> >		("gsftrkscharge"	).setBranchAlias("gsftrks_charge"     );	// charge						


  
  produces<vector<int> >		  ("gsftrksnlayers"	).setBranchAlias("gsftrks_nlayers"    );
  produces<vector<int> >		  ("gsftrksnlayers3D"	).setBranchAlias("gsftrks_nlayers3D"  );
  produces<vector<int> >		  ("gsftrksnlayersLost"	).setBranchAlias("gsftrks_nlayersLost");

   //Hit Pattern information
  produces<vector<LorentzVector> >	  ("gsftrksinnerposition"     ).setBranchAlias("gsftrks_inner_position"         );
  produces<vector<LorentzVector> >	  ("gsftrksouterposition"     ).setBranchAlias("gsftrks_outer_position"         );
  produces<vector<int> >		  ("gsftrksvalidpixelhits"    ).setBranchAlias("gsftrks_valid_pixelhits"        );
  produces<vector<int> >		  ("gsftrkslostpixelhits"     ).setBranchAlias("gsftrks_lost_pixelhits"         );
  produces<vector<int> >		  ("gsftrkslayer1sizerphi"    ).setBranchAlias("gsftrks_layer1_sizerphi"        ); 
  produces<vector<int> >		  ("gsftrkslayer1sizerz"      ).setBranchAlias("gsftrks_layer1_sizerz"          ); 
  produces<vector<float> >		  ("gsftrkslayer1charge"      ).setBranchAlias("gsftrks_layer1_charge"          ); 
  produces<vector<int> >		  ("gsftrkslayer1det"         ).setBranchAlias("gsftrks_layer1_det"             );
  produces<vector<int> >		  ("gsftrkslayer1layer"       ).setBranchAlias("gsftrks_layer1_layer"           ); 
  produces<vector<int> >		  ("gsftrksexpinnerlayers"    ).setBranchAlias("gsftrks_exp_innerlayers"        );
  produces<vector<int> >		  ("gsftrksexpouterlayers"    ).setBranchAlias("gsftrks_exp_outerlayers"        );   

  gsftracksInputTag_ = iConfig.getParameter<edm::InputTag>("gsftracksInputTag");
  beamSpotTag_    = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");


}

void GSFTrackMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  std::auto_ptr<vector<LorentzVector> >	gsftrks_p4		(new vector<LorentzVector>	);
  std::auto_ptr<vector<LorentzVector> >	gsftrks_vertex_p4	(new vector<LorentzVector>	);
  std::auto_ptr<vector<LorentzVector> > gsftrks_outer_p4		(new vector<LorentzVector>      );
  std::auto_ptr<vector<float> >		gsftrks_d0		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_d0corr		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_d0corrPhi	(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_z0		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_z0corr		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_chi2		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_ndof		(new vector<float>		);      
  std::auto_ptr<vector<int> >		gsftrks_validHits	(new vector<int>		);        
  std::auto_ptr<vector<int> >		gsftrks_lostHits		(new vector<int>		);        
  std::auto_ptr<vector<float> >		gsftrks_d0Err		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_z0Err		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_ptErr		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_etaErr		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_phiErr		(new vector<float>		);      
  std::auto_ptr<vector<float> >		gsftrks_d0phiCov		(new vector<float>		);      
  std::auto_ptr<vector<int> >		gsftrks_charge		(new vector<int>		);        


   //HitPattern information
  //
  std::auto_ptr<vector<LorentzVector> >gsftrks_inner_position		(new vector<LorentzVector>	);
  std::auto_ptr<vector<LorentzVector> >gsftrks_outer_position		(new vector<LorentzVector>	);
  std::auto_ptr<vector<int> >	       gsftrks_valid_pixelhits		(new vector<int>	        ); 
  std::auto_ptr<vector<int> >	       gsftrks_lost_pixelhits		(new vector<int>        	); 
  std::auto_ptr<vector<int> >	       gsftrks_layer1_sizerphi		(new vector<int>	        ); 
  std::auto_ptr<vector<int> >	       gsftrks_layer1_sizerz		(new vector<int>	        ); 
  std::auto_ptr<vector<float> >	       gsftrks_layer1_charge		(new vector<float>	        );
  std::auto_ptr<vector<int> >	       gsftrks_layer1_det		(new vector<int>	        );
  std::auto_ptr<vector<int> >	       gsftrks_layer1_layer		(new vector<int>		);
  std::auto_ptr<vector<int> >	       gsftrks_exp_innerlayers		(new vector<int>		); 
  std::auto_ptr<vector<int> >	       gsftrks_exp_outerlayers		(new vector<int>		); 

  std::auto_ptr<vector<int> > gsftrks_nlayers			(new vector<int>		);
  std::auto_ptr<vector<int> > gsftrks_nlayers3D			(new vector<int>		);
  std::auto_ptr<vector<int> > gsftrks_nlayersLost		(new vector<int>		);

  // get tracks
  Handle<edm::View<reco::GsfTrack> > track_h;
  iEvent.getByLabel(gsftracksInputTag_, track_h);

  if( !track_h.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve track collection";
    edm::LogInfo("OutputInfo") << " GSFTrackMaker cannot continue...!";
    return;
  }

  edm::ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);

  if( !theMagField.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve the magnetic field";
    edm::LogInfo("OutputInfo") << " GSFTrackMaker cannot continue...!";
    return;
  }

  const MagneticField* bf = theMagField.product();
  
  edm::Handle<LorentzVector> beamSpotH;
  iEvent.getByLabel(beamSpotTag_, beamSpotH);
  
  const Point beamSpot = beamSpotH.isValid() ? Point(beamSpotH->x(), beamSpotH->y(), beamSpotH->z()) : Point(0, 0, 0);
  
  //get tracker geometry   
                                                                                                                                                                  
  edm::ESHandle<TrackerGeometry> theG;
  iSetup.get<TrackerDigiGeometryRecord>().get(theG);
  

  edm::View<reco::GsfTrack>::const_iterator tracks_end = track_h->end();
  
  for (edm::View<reco::GsfTrack>::const_iterator i = track_h->begin(); i != tracks_end; ++i) {

    gsftrks_p4           ->push_back( LorentzVector( i->px(), i->py(), i->pz(), i->p() )       );
    gsftrks_vertex_p4    ->push_back( LorentzVector(i->vx(),i->vy(), i->vz(), 0.)              );
    gsftrks_d0           ->push_back( i->d0()                                                  );
    gsftrks_z0           ->push_back( i->dz()                                                  );
    						           
    double corrd0 = beamSpotH.isValid() ? -1 * ( i->dxy(beamSpot) ) : i->d0();		           
    gsftrks_d0corr       ->push_back( corrd0                                                   );
    						           
    double corrd0phi = atan2( -1 * corrd0 * sin( i->phi() ), corrd0 * cos( i->phi() ) );           
    gsftrks_d0corrPhi    ->push_back( corrd0phi                                                );
    						           
    double corrz0 = beamSpotH.isValid() ? i->dz(beamSpot) : i->dz();			           
    gsftrks_z0corr       ->push_back( corrz0                                                   );
    						           
    gsftrks_chi2         ->push_back( i->chi2()                                                );
    gsftrks_ndof         ->push_back( i->ndof()                                                );
    gsftrks_validHits    ->push_back( i->numberOfValidHits()                                   );
    gsftrks_lostHits     ->push_back( i->numberOfLostHits()                                    );
    gsftrks_d0Err        ->push_back( i->d0Error()                                             );
    gsftrks_z0Err        ->push_back( i->dzError()                                             );
    gsftrks_ptErr        ->push_back( i->ptError()                                             );
    gsftrks_etaErr       ->push_back( i->etaError()                                            );
    gsftrks_phiErr       ->push_back( i->phiError()                                            );
    gsftrks_d0phiCov     ->push_back( -i->covariance(TrackBase::i_phi, TrackBase::i_dxy)       );  // minus sign because we want cov(d0, phi) = cov(-dxy, phi)
    gsftrks_charge       ->push_back( i->charge()                                              );

    

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
      gsftrks_outer_p4->push_back( LorentzVector( tsos.globalMomentum().x(),
						      tsos.globalMomentum().y(),
						      tsos.globalMomentum().z(),
						      tsos.globalMomentum().mag() ) );
    }
    else {
      gsftrks_outer_p4->push_back( LorentzVector( 999., 0., 22004439., 22004440.) );
   }

    /////hit pattern
    gsftrks_inner_position ->push_back(LorentzVector(i->innerPosition().x(), i->innerPosition().y() , i->innerPosition().z(), 0 ));
    gsftrks_outer_position ->push_back(LorentzVector(i->outerPosition().x(), i->outerPosition().y() , i->outerPosition().z(), 0 ));
    const reco::HitPattern& pattern = i->hitPattern();
    const reco::HitPattern& p_inner = i->trackerExpectedHitsInner();
    const reco::HitPattern& p_outer = i->trackerExpectedHitsOuter();
    gsftrks_exp_innerlayers    -> push_back(p_inner.numberOfHits());
    gsftrks_exp_outerlayers    -> push_back(p_outer.numberOfHits());
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


    for(trackingRecHit_iterator ihit = i->recHitsBegin(); 
	ihit != i->recHitsEnd(); ++ihit){
      if(i_layer > 1) break;
      int k = ihit-i->recHitsBegin();
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
	  gsftrks_layer1_sizerphi ->push_back(pixel_sizeX);
	  gsftrks_layer1_sizerz   ->push_back(pixel_sizeY);
	  gsftrks_layer1_charge   ->push_back(pixel_charge);
	  gsftrks_layer1_det      ->push_back(det);
	  gsftrks_layer1_layer    ->push_back(layer);
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
	      gsftrks_layer1_sizerphi ->push_back(cluster_size);
	      gsftrks_layer1_sizerz   ->push_back(0);
	    }

	  else
	    {
	      gsftrks_layer1_sizerphi ->push_back(0);
	      gsftrks_layer1_sizerz   ->push_back(cluster_size);
	    } 
	  gsftrks_layer1_charge   ->push_back(cluster_charge);
	  gsftrks_layer1_det      ->push_back(det);
	  gsftrks_layer1_layer    ->push_back(layer);
	  i_layer++;
	}
      }
    }
    gsftrks_valid_pixelhits ->push_back(pattern.numberOfValidPixelHits());
    gsftrks_lost_pixelhits ->push_back(pattern.numberOfLostPixelHits());
    
    
    // *****************************************************
     gsftrks_nlayers    ->push_back( i->hitPattern().trackerLayersWithMeasurement() );
     gsftrks_nlayers3D  ->push_back( i->hitPattern().pixelLayersWithMeasurement() + i->hitPattern().numberOfValidStripLayersWithMonoAndStereo() );
     gsftrks_nlayersLost->push_back( i->hitPattern().trackerLayersWithoutMeasurement() );

  }

  // store vectors
  iEvent.put(gsftrks_p4			, "gsftrksp4"			);
  iEvent.put(gsftrks_vertex_p4		, "gsftrksvertexp4"		);
  iEvent.put(gsftrks_outer_p4		, "gsftrksouterp4"		);
  iEvent.put(gsftrks_d0			, "gsftrksd0"			);
  iEvent.put(gsftrks_d0corr		, "gsftrksd0corr"		);
  iEvent.put(gsftrks_d0corrPhi		, "gsftrksd0corrPhi"		);
  iEvent.put(gsftrks_z0			, "gsftrksz0"			);
  iEvent.put(gsftrks_z0corr		, "gsftrksz0corr"		);
  iEvent.put(gsftrks_chi2		, "gsftrkschi2"			);
  iEvent.put(gsftrks_ndof		, "gsftrksndof"			);
  iEvent.put(gsftrks_validHits		, "gsftrksvalidHits"		);
  iEvent.put(gsftrks_lostHits		, "gsftrkslostHits"		);
  iEvent.put(gsftrks_d0Err		, "gsftrksd0Err"		);
  iEvent.put(gsftrks_z0Err		, "gsftrksz0Err"		);
  iEvent.put(gsftrks_ptErr		, "gsftrksptErr"		);
  iEvent.put(gsftrks_etaErr		, "gsftrksetaErr"		);
  iEvent.put(gsftrks_phiErr		, "gsftrksphiErr"		);
  iEvent.put(gsftrks_d0phiCov		, "gsftrksd0phiCov"		);
  iEvent.put(gsftrks_charge		, "gsftrkscharge"		);
  iEvent.put(gsftrks_nlayers		, "gsftrksnlayers"		);
  iEvent.put(gsftrks_nlayers3D		, "gsftrksnlayers3D"		);
  iEvent.put(gsftrks_nlayersLost	, "gsftrksnlayersLost"		);

  //Hit Pattern Information
  iEvent.put(gsftrks_inner_position	, "gsftrksinnerposition"	);
  iEvent.put(gsftrks_outer_position	, "gsftrksouterposition"	);
  iEvent.put(gsftrks_valid_pixelhits	, "gsftrksvalidpixelhits"	);
  iEvent.put(gsftrks_lost_pixelhits	, "gsftrkslostpixelhits"	);
  iEvent.put(gsftrks_layer1_layer	, "gsftrkslayer1layer"		);
  iEvent.put(gsftrks_layer1_sizerphi	, "gsftrkslayer1sizerphi"	);
  iEvent.put(gsftrks_layer1_sizerz	, "gsftrkslayer1sizerz"		);
  iEvent.put(gsftrks_layer1_charge	, "gsftrkslayer1charge"		);
  iEvent.put(gsftrks_layer1_det		, "gsftrkslayer1det"		);
  iEvent.put(gsftrks_exp_innerlayers	, "gsftrksexpinnerlayers"	);
  iEvent.put(gsftrks_exp_outerlayers	, "gsftrksexpouterlayers"	);
  
}

// ------------ method called once each job just before starting event loop  ------------
void GSFTrackMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GSFTrackMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(GSFTrackMaker);
