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
// $Id: TrackMaker.cc,v 1.4 2008/07/22 06:08:21 fgolf Exp $
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
#include "CMS2/NtupleMaker/interface/TrackMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
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
using std::vector;
//
// class decleration
//

//
// constructors and destructor
//
TrackMaker::TrackMaker(const edm::ParameterSet& iConfig)
{
     // stream mu track quantities
     produces<vector<LorentzVector> >	("trkstrkp4"		).setBranchAlias("trks_trk_p4"       	);	// track p4						
     produces<vector<float> >		("trksd0"		).setBranchAlias("trks_d0"           	);	// impact parameter at the point of closest approach	
     produces<vector<float> >		("trksz0"		).setBranchAlias("trks_z0"           	);	// z position of the point of closest approach		
     produces<vector<float> >		("trksvertexphi"	).setBranchAlias("trks_vertexphi"    	);	// phi angle of the point of closest approach		
     produces<vector<float> >		("trkschi2"		).setBranchAlias("trks_chi2"         	);	// chi2 of the silicon tracker fit			
     produces<vector<float> >		("trksndof"		).setBranchAlias("trks_ndof"         	);	// number of degrees of freedom of the fit		
     produces<vector<int> >		("trksvalidHits"	).setBranchAlias("trks_validHits"    	);	// number of used hits in the fit			
     produces<vector<int> >		("trkslostHits"		).setBranchAlias("trks_lostHits"     	);	// number of lost hits in the fit			
     produces<vector<float> >		("trksd0Err"		).setBranchAlias("trks_d0Err"        	);	// error on the impact parameter			
     produces<vector<float> >		("trksz0Err"		).setBranchAlias("trks_z0Err"        	);	// error on z position of the point of closest approach	
     produces<vector<float> >		("trksptErr"		).setBranchAlias("trks_ptErr"        	);	// track Pt error					
     produces<vector<float> >		("trksetaErr"		).setBranchAlias("trks_etaErr"       	);	// track eta error					
     produces<vector<float> >		("trksphiErr"		).setBranchAlias("trks_phiErr"       	);	// track phi error					
     produces<vector<int> >		("trkscharge"		).setBranchAlias("trks_charge"       	);	// charge						
     produces<vector<float> >		("trksouterPhi"		).setBranchAlias("trks_outerPhi"     	);	// phi angle of the outermost point in tracker		
     produces<vector<float> >		("trksouterEta"		).setBranchAlias("trks_outerEta"     	);	// eta angle of the outermost point in tracker		

     tracksInputTag = iConfig.getParameter<edm::InputTag>("tracksInputTag");
}

void TrackMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<LorentzVector> >	vector_trks_p4		(new vector<LorentzVector>	);
     std::auto_ptr<vector<LorentzVector> >	vector_trks_trk_p4	(new vector<LorentzVector>	);
     std::auto_ptr<vector<float> >		vector_trks_d0		(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_trks_z0		(new vector<float>		);      
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
     // get tracks
     Handle<edm::View<reco::Track> > track_h;
     iEvent.getByLabel(tracksInputTag, track_h);      // change this in the future
     edm::View<reco::Track>::const_iterator tracks_end = track_h->end();
     // get magnetic field
     edm::ESHandle<MagneticField> theMagField;
     iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
     const MagneticField* bf = theMagField.product();

     for (edm::View<reco::Track>::const_iterator i = track_h->begin(); 
	  i != tracks_end; ++i) {
	  // fill vectors
          vector_trks_trk_p4       ->push_back(	LorentzVector( i->px(), i->py(), i->pz(), i->p() ) );
	  vector_trks_d0           ->push_back( i->d0() );
	  vector_trks_z0           ->push_back(	i->dz() );
	  vector_trks_vertexphi    ->push_back( atan2( i->vy(), i->vx() ) );
	  vector_trks_chi2         ->push_back(	i->chi2() );
	  vector_trks_ndof         ->push_back(	i->ndof() );
	  vector_trks_validHits    ->push_back(	i->numberOfValidHits() );
	  vector_trks_lostHits     ->push_back(	i->numberOfLostHits() );
	  vector_trks_d0Err        ->push_back(	i->d0Error() );
	  vector_trks_z0Err        ->push_back(	i->dzError() );
	  vector_trks_ptErr        ->push_back(	i->ptError() );
	  vector_trks_etaErr       ->push_back(	i->etaError() );
	  vector_trks_phiErr       ->push_back(	i->phiError() );
	  vector_trks_charge       ->push_back(	i->charge() );
	  
	  GlobalPoint  tpVertex ( i->vx(), i->vy(), i->vz() );
	  GlobalVector tpMomentum ( i->px(), i->py(), i->pz() );
	  int tpCharge ( i->charge() );

	  FreeTrajectoryState fts ( tpVertex, tpMomentum, tpCharge, bf);

	  const float zdist = 314.;

	  const float radius = 130.;

	  const float corner = 1.479;

	  Plane::PlanePointer lendcap = Plane::build( Plane::PositionType (0, 0, -zdist), Plane::RotationType () );
	  Plane::PlanePointer rendcap = Plane::build( Plane::PositionType (0, 0, zdist), Plane::RotationType () );

	  Cylinder::CylinderPointer barrel = Cylinder::build( Cylinder::PositionType (0, 0, 0), Cylinder::RotationType (), radius);

	  AnalyticalPropagator myAP (bf, alongMomentum, 2*M_PI);

	  TrajectoryStateOnSurface tsos;

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
	    vector_trks_outerPhi->push_back( tsos.globalPosition().phi() );
	    vector_trks_outerEta->push_back( tsos.globalPosition().eta() );
	  }
	  else {
	    vector_trks_outerPhi->push_back( -999. );
	    vector_trks_outerEta->push_back( -999. );
	  }
     }
     // store vectors
     iEvent.put(vector_trks_trk_p4       , "trkstrkp4"             );
     iEvent.put(vector_trks_d0           , "trksd0"                );
     iEvent.put(vector_trks_z0           , "trksz0"                );
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
     iEvent.put(vector_trks_outerPhi     , "trksouterPhi"          );
     iEvent.put(vector_trks_outerEta     , "trksouterEta"          );
}

// ------------ method called once each job just before starting event loop  ------------
void 
TrackMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackMaker);
