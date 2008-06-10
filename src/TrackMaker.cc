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
// $Id: TrackMaker.cc,v 1.2 2008/06/10 20:06:15 jmuelmen Exp $
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
     produces<vector<LorentzVector> >	("musp4"		).setBranchAlias("mus_p4"           	);	// candidate p4						
     produces<vector<LorentzVector> >	("mustrkp4"		).setBranchAlias("mus_trk_p4"       	);	// track p4						
     produces<vector<float> >		("musd0"		).setBranchAlias("mus_d0"           	);	// impact parameter at the point of closest approach	
     produces<vector<float> >		("musz0"		).setBranchAlias("mus_z0"           	);	// z position of the point of closest approach		
     produces<vector<float> >		("musvertexphi"		).setBranchAlias("mus_vertexphi"    	);	// phi angle of the point of closest approach		
     produces<vector<float> >		("muschi2"		).setBranchAlias("mus_chi2"         	);	// chi2 of the silicon tracker fit			
     produces<vector<float> >		("musndof"		).setBranchAlias("mus_ndof"         	);	// number of degrees of freedom of the fit		
     produces<vector<int> >		("musvalidHits"		).setBranchAlias("mus_validHits"    	);	// number of used hits in the fit			
     produces<vector<int> >		("muslostHits"		).setBranchAlias("mus_lostHits"     	);	// number of lost hits in the fit			
     produces<vector<float> >		("musd0Err"		).setBranchAlias("mus_d0Err"        	);	// error on the impact parameter			
     produces<vector<float> >		("musz0Err"		).setBranchAlias("mus_z0Err"        	);	// error on z position of the point of closest approach	
     produces<vector<float> >		("musptErr"		).setBranchAlias("mus_ptErr"        	);	// track Pt error					
     produces<vector<float> >		("musetaErr"		).setBranchAlias("mus_etaErr"       	);	// track eta error					
     produces<vector<float> >		("musphiErr"		).setBranchAlias("mus_phiErr"       	);	// track phi error					
     produces<vector<LorentzVector> >	("musmcp4"		).setBranchAlias("mus_mc_p4"        	);	// p4 of matched MC particle				
     produces<vector<int> >		("musmcid"		).setBranchAlias("mus_mc_id"        	);	// PDG id of matched MC particle			
     produces<vector<int> >		("muscharge"		).setBranchAlias("mus_charge"       	);	// charge						
     produces<vector<int> >		("musmcmotherid"	).setBranchAlias("mus_mc_motherid"  	);	// PDG id of the mother of the particle			
     produces<vector<float> >		("musouterPhi"		).setBranchAlias("mus_outerPhi"     	);	// phi angle of the outermost point in tracker		
     produces<vector<float> >		("musouterEta"		).setBranchAlias("mus_outerEta"     	);	// eta angle of the outermost point in tracker		
}

void TrackMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<LorentzVector> >	vector_mus_p4		(new vector<LorentzVector>	);
     std::auto_ptr<vector<LorentzVector> >	vector_mus_trk_p4	(new vector<LorentzVector>	);
     std::auto_ptr<vector<float> >		vector_mus_d0		(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_mus_z0		(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_mus_vertexphi	(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_mus_chi2		(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_mus_ndof		(new vector<float>		);      
     std::auto_ptr<vector<int> >		vector_mus_validHits	(new vector<int>		);        
     std::auto_ptr<vector<int> >		vector_mus_lostHits	(new vector<int>		);        
     std::auto_ptr<vector<float> >		vector_mus_d0Err	(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_mus_z0Err	(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_mus_ptErr	(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_mus_etaErr	(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_mus_phiErr	(new vector<float>		);      
     std::auto_ptr<vector<LorentzVector> >	vector_mus_mc_p4	(new vector<LorentzVector>	);
     std::auto_ptr<vector<int> >		vector_mus_mc_id	(new vector<int>		);        
     std::auto_ptr<vector<int> >		vector_mus_charge	(new vector<int>		);        
     std::auto_ptr<vector<int> >		vector_mus_mc_motherid	(new vector<int>		);        
     std::auto_ptr<vector<float> >		vector_mus_outerPhi	(new vector<float>		);      
     std::auto_ptr<vector<float> >		vector_mus_outerEta	(new vector<float>		);      
     // get tracks
     Handle<edm::View<reco::Track> > track_h;
     iEvent.getByLabel("ctfWithMaterialTracks", track_h);      // change this in the future
     edm::View<reco::Track>::const_iterator tracks_end = track_h->end();
     for (edm::View<reco::Track>::const_iterator i = track_h->begin(); 
	  i != tracks_end; ++i) {
	  // fill vectors
#if 0
	  vector_mus_p4           ->push_back(	);
	  vector_mus_trk_p4       ->push_back(	);
	  vector_mus_d0           ->push_back(	);
	  vector_mus_z0           ->push_back(	);
	  vector_mus_vertexphi    ->push_back(	);
	  vector_mus_chi2         ->push_back(	);
	  vector_mus_ndof         ->push_back(	);
	  vector_mus_validHits    ->push_back(	);
	  vector_mus_lostHits     ->push_back(	);
	  vector_mus_d0Err        ->push_back(	);
	  vector_mus_z0Err        ->push_back(	);
	  vector_mus_ptErr        ->push_back(	);
	  vector_mus_etaErr       ->push_back(	);
	  vector_mus_phiErr       ->push_back(	);
	  vector_mus_mc_p4        ->push_back(	);
	  vector_mus_mc_id        ->push_back(	);
	  vector_mus_charge       ->push_back(	);
	  vector_mus_mc_motherid  ->push_back(	);
	  vector_mus_outerPhi     ->push_back(	);
	  vector_mus_outerEta     ->push_back(	);
#endif
     }
     // store vectors
     iEvent.put(vector_mus_p4           , "musp4"                );
     iEvent.put(vector_mus_trk_p4       , "mustrkp4"             );
     iEvent.put(vector_mus_d0           , "musd0"                );
     iEvent.put(vector_mus_z0           , "musz0"                );
     iEvent.put(vector_mus_vertexphi    , "musvertexphi"         );
     iEvent.put(vector_mus_chi2         , "muschi2"              );
     iEvent.put(vector_mus_ndof         , "musndof"              );
     iEvent.put(vector_mus_validHits    , "musvalidHits"         );
     iEvent.put(vector_mus_lostHits     , "muslostHits"          );
     iEvent.put(vector_mus_d0Err        , "musd0Err"             );
     iEvent.put(vector_mus_z0Err        , "musz0Err"             );
     iEvent.put(vector_mus_ptErr        , "musptErr"             );
     iEvent.put(vector_mus_etaErr       , "musetaErr"            );
     iEvent.put(vector_mus_phiErr       , "musphiErr"            );
     iEvent.put(vector_mus_mc_p4        , "musmcp4"              );
     iEvent.put(vector_mus_mc_id        , "musmcid"              );
     iEvent.put(vector_mus_charge       , "muscharge"            );
     iEvent.put(vector_mus_mc_motherid  , "musmcmotherid"        );
     iEvent.put(vector_mus_outerPhi     , "musouterPhi"          );
     iEvent.put(vector_mus_outerEta     , "musouterEta"          );
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
