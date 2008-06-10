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
// $Id: MuonMaker.cc,v 1.2 2008/06/10 19:42:37 jmuelmen Exp $
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
#include "RecoMuon/MuonIdentification/interface/IdGlobalFunctions.h"

typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

MuonMaker::MuonMaker(const edm::ParameterSet& iConfig)
{
     // mu track quantities
     produces<vector<LorentzVector> >	("musp4"		).setBranchAlias("mus_p4"           	);	// candidate p4						
     produces<vector<LorentzVector> >	("mustrk_p4"		).setBranchAlias("mus_trk_p4"       	);	// track p4						
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
     produces<vector<LorentzVector> >	("musmc_p4"		).setBranchAlias("mus_mc_p4"        	);	// p4 of matched MC particle				
     produces<vector<int> >		("musmc_id"		).setBranchAlias("mus_mc_id"        	);	// PDG id of matched MC particle			
     produces<vector<int> >		("muscharge"		).setBranchAlias("mus_charge"       	);	// charge						
     produces<vector<int> >		("musmc_motherid"	).setBranchAlias("mus_mc_motherid"  	);	// PDG id of the mother of the particle			
     produces<vector<float> >		("musouterPhi"		).setBranchAlias("mus_outerPhi"     	);	// phi angle of the outermost point in tracker		
     produces<vector<float> >		("musouterEta"		).setBranchAlias("mus_outerEta"     	);	// eta angle of the outermost point in tracker		
     // muon quantities
     produces<vector<int> >	("musnmatches"			).setBranchAlias("mus_nmatches"                     	);
     produces<vector<float> >	("museem"			).setBranchAlias("mus_e_em"                         	);
     produces<vector<float> >	("musehad"			).setBranchAlias("mus_e_had"                        	);
     produces<vector<float> >	("museho"			).setBranchAlias("mus_e_ho"                         	);
     produces<vector<float> >	("museemS9"			).setBranchAlias("mus_e_emS9"                       	);
     produces<vector<float> >	("musehadS9"			).setBranchAlias("mus_e_hadS9"                      	);
     produces<vector<float> >	("musehoS9"			).setBranchAlias("mus_e_hoS9"                       	);
     produces<vector<float> >	("musiso03sumPt"		).setBranchAlias("mus_iso03_sumPt"                  	);
     produces<vector<float> >	("musiso03emEt"			).setBranchAlias("mus_iso03_emEt"                   	);
     produces<vector<float> >	("musiso03hadEt"		).setBranchAlias("mus_iso03_hadEt"                  	);
     produces<vector<float> >	("musiso03hoEt"			).setBranchAlias("mus_iso03_hoEt"                   	);
     produces<vector<int> >	("musiso03ntrk"			).setBranchAlias("mus_iso03_ntrk"                   	);
     produces<vector<float> >	("musiso05sumPt"		).setBranchAlias("mus_iso05_sumPt"                  	);
     produces<vector<float> >	("musiso05emEt"			).setBranchAlias("mus_iso05_emEt"                   	);
     produces<vector<float> >	("musiso05hadEt"		).setBranchAlias("mus_iso05_hadEt"                  	);
     produces<vector<float> >	("musiso05hoEt"			).setBranchAlias("mus_iso05_hoEt"                   	);
     produces<vector<int> >	("musiso05ntrk"			).setBranchAlias("mus_iso05_ntrk"                   	);
     produces<vector<int> >	("muspidTMLastStationLoose"	).setBranchAlias("mus_pid_TMLastStationLoose"       	);
     produces<vector<int> >	("muspidTMLastStationTight"	).setBranchAlias("mus_pid_TMLastStationTight"       	);
     produces<vector<int> >	("muspidTM2DCompatibilityLoose"	).setBranchAlias("mus_pid_TM2DCompatibilityLoose"   	);
     produces<vector<int> >	("muspidTM2DCompatibilityTight"	).setBranchAlias("mus_pid_TM2DCompatibilityTight"   	);
     produces<vector<float> >	("musgfitchi2"			).setBranchAlias("mus_gfit_chi2"                    	);
     produces<vector<float> >	("musgfitndof"			).setBranchAlias("mus_gfit_ndof"                    	);
     produces<vector<int> >	("musgfitvalidHits"		).setBranchAlias("mus_gfit_validHits"               	);
}

void MuonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
     std::auto_ptr<vector<int> >		vector_mus_nmatches			(new vector<int>  	);
     std::auto_ptr<vector<float> >		vector_mus_e_em				(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_e_had			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_e_ho				(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_e_emS9			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_e_hadS9			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_e_hoS9			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_iso03_sumPt			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_iso03_emEt			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_iso03_hadEt			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_iso03_hoEt			(new vector<float>	);
     std::auto_ptr<vector<int> >		vector_mus_iso03_ntrk			(new vector<int>  	);
     std::auto_ptr<vector<float> >		vector_mus_iso05_sumPt			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_iso05_emEt			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_iso05_hadEt			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_iso05_hoEt			(new vector<float>	);
     std::auto_ptr<vector<int> >		vector_mus_iso05_ntrk			(new vector<int>  	);
     std::auto_ptr<vector<int> >		vector_mus_pid_TMLastStationLoose	(new vector<int>  	);
     std::auto_ptr<vector<int> >		vector_mus_pid_TMLastStationTight	(new vector<int>  	);
     std::auto_ptr<vector<int> >		vector_mus_pid_TM2DCompatibilityLoose	(new vector<int>  	);
     std::auto_ptr<vector<int> >		vector_mus_pid_TM2DCompatibilityTight	(new vector<int>  	);
     std::auto_ptr<vector<float> >		vector_mus_gfit_chi2			(new vector<float>	);
     std::auto_ptr<vector<float> >		vector_mus_gfit_ndof			(new vector<float>	);
     std::auto_ptr<vector<int> >		vector_mus_gfit_validHits		(new vector<int>  	);
     // get muons
     Handle<edm::View<reco::Muon> > muon_h;
     iEvent.getByLabel("allLayer1TopMuons", muon_h);      // change this in the future
     edm::View<reco::Muon>::const_iterator muons_end = muon_h->end();
     for (edm::View<reco::Muon>::const_iterator muon = muon_h->begin(); 
	  muon != muons_end; ++muon) {
	  // fill vectors
#if 0
	  // track information will be copied from TrackMaker later...
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
	  vector_mus_nmatches                     ->push_back(muon->isMatchesValid() ? muon->numberOfMatches()			: -999	);
	  vector_mus_e_em                         ->push_back(muon->isEnergyValid() ? muon->getCalEnergy().em			: -999	);
	  vector_mus_e_had                        ->push_back(muon->isEnergyValid() ? muon->getCalEnergy().had			: -999	);
	  vector_mus_e_ho                         ->push_back(muon->isEnergyValid() ? muon->getCalEnergy().ho			: -999	);
	  vector_mus_e_emS9                       ->push_back(muon->isEnergyValid() ? muon->getCalEnergy().emS9			: -999	);
	  vector_mus_e_hadS9                      ->push_back(muon->isEnergyValid() ? muon->getCalEnergy().hadS9		: -999	);
	  vector_mus_e_hoS9                       ->push_back(muon->isEnergyValid() ? muon->getCalEnergy().hoS9			: -999	);
	  vector_mus_iso03_sumPt                  ->push_back(muon->isIsolationValid() ? muon->getIsolationR03().sumPt		: -999	);
	  vector_mus_iso03_emEt                   ->push_back(muon->isIsolationValid() ? muon->getIsolationR03().emEt		: -999	);
	  vector_mus_iso03_hadEt                  ->push_back(muon->isIsolationValid() ? muon->getIsolationR03().hadEt		: -999	);
	  vector_mus_iso03_hoEt                   ->push_back(muon->isIsolationValid() ? muon->getIsolationR03().hoEt		: -999	);
	  vector_mus_iso03_ntrk                   ->push_back(muon->isIsolationValid() ? muon->getIsolationR03().nTracks	: -999	);
	  vector_mus_iso05_sumPt                  ->push_back(muon->isIsolationValid() ? muon->getIsolationR05().sumPt		: -999	);
	  vector_mus_iso05_emEt                   ->push_back(muon->isIsolationValid() ? muon->getIsolationR05().emEt		: -999	);
	  vector_mus_iso05_hadEt                  ->push_back(muon->isIsolationValid() ? muon->getIsolationR05().hadEt		: -999	);
	  vector_mus_iso05_hoEt                   ->push_back(muon->isIsolationValid() ? muon->getIsolationR05().hoEt		: -999	);
	  vector_mus_iso05_ntrk                   ->push_back(muon->isIsolationValid() ? muon->getIsolationR05().nTracks	: -999	);
	  vector_mus_pid_TMLastStationLoose       ->push_back(muon->isMatchesValid() ? muonid::isGoodMuon(*muon,muonid::TMLastStationLoose)     : -999	);
	  vector_mus_pid_TMLastStationTight       ->push_back(muon->isMatchesValid() ? muonid::isGoodMuon(*muon,muonid::TMLastStationTight)     : -999	);
	  vector_mus_pid_TM2DCompatibilityLoose   ->push_back(muon->isMatchesValid() ? muonid::isGoodMuon(*muon,muonid::TM2DCompatibilityLoose)	: -999	);
	  vector_mus_pid_TM2DCompatibilityTight   ->push_back(muon->isMatchesValid() ? muonid::isGoodMuon(*muon,muonid::TM2DCompatibilityTight)	: -999	);
	  vector_mus_gfit_chi2                    ->push_back(muon->combinedMuon().isNonnull() ? muon->combinedMuon()->chi2() 			: -999	);
	  vector_mus_gfit_ndof                    ->push_back(muon->combinedMuon().isNonnull() ? muon->combinedMuon()->ndof() 			: -999	);
     	  vector_mus_gfit_validHits               ->push_back(muon->combinedMuon().isNonnull() ? muon->combinedMuon()->numberOfValidHits() 	: -999	);
     }
     // store vectors
     iEvent.put(vector_mus_p4           , "mus_p4"               );
     iEvent.put(vector_mus_trk_p4       , "mus_trk_p4"           );
     iEvent.put(vector_mus_d0           , "mus_d0"               );
     iEvent.put(vector_mus_z0           , "mus_z0"               );
     iEvent.put(vector_mus_vertexphi    , "mus_vertexphi"        );
     iEvent.put(vector_mus_chi2         , "mus_chi2"             );
     iEvent.put(vector_mus_ndof         , "mus_ndof"             );
     iEvent.put(vector_mus_validHits    , "mus_validHits"        );
     iEvent.put(vector_mus_lostHits     , "mus_lostHits"         );
     iEvent.put(vector_mus_d0Err        , "mus_d0Err"            );
     iEvent.put(vector_mus_z0Err        , "mus_z0Err"            );
     iEvent.put(vector_mus_ptErr        , "mus_ptErr"            );
     iEvent.put(vector_mus_etaErr       , "mus_etaErr"           );
     iEvent.put(vector_mus_phiErr       , "mus_phiErr"           );
     iEvent.put(vector_mus_mc_p4        , "mus_mc_p4"            );
     iEvent.put(vector_mus_mc_id        , "mus_mc_id"            );
     iEvent.put(vector_mus_charge       , "mus_charge"           );
     iEvent.put(vector_mus_mc_motherid  , "mus_mc_motherid"      );
     iEvent.put(vector_mus_outerPhi     , "mus_outerPhi"         );
     iEvent.put(vector_mus_outerEta     , "mus_outerEta"         );
     iEvent.put(vector_mus_nmatches			, "musnmatches"                  );
     iEvent.put(vector_mus_e_em				, "museem"                       );
     iEvent.put(vector_mus_e_had			, "musehad"                      );
     iEvent.put(vector_mus_e_ho				, "museho"                       );
     iEvent.put(vector_mus_e_emS9			, "museemS9"                     );
     iEvent.put(vector_mus_e_hadS9			, "musehadS9"                    );
     iEvent.put(vector_mus_e_hoS9			, "musehoS9"                     );
     iEvent.put(vector_mus_iso03_sumPt			, "musiso03sumPt"                );
     iEvent.put(vector_mus_iso03_emEt			, "musiso03emEt"                 );
     iEvent.put(vector_mus_iso03_hadEt			, "musiso03hadEt"                );
     iEvent.put(vector_mus_iso03_hoEt			, "musiso03hoEt"                 );
     iEvent.put(vector_mus_iso03_ntrk			, "musiso03ntrk"                 );
     iEvent.put(vector_mus_iso05_sumPt			, "musiso05sumPt"                );
     iEvent.put(vector_mus_iso05_emEt			, "musiso05emEt"                 );
     iEvent.put(vector_mus_iso05_hadEt			, "musiso05hadEt"                );
     iEvent.put(vector_mus_iso05_hoEt			, "musiso05hoEt"                 );
     iEvent.put(vector_mus_iso05_ntrk			, "musiso05ntrk"                 );
     iEvent.put(vector_mus_pid_TMLastStationLoose	, "muspidTMLastStationLoose"     );
     iEvent.put(vector_mus_pid_TMLastStationTight	, "muspidTMLastStationTight"     );
     iEvent.put(vector_mus_pid_TM2DCompatibilityLoose	, "muspidTM2DCompatibilityLoose" );
     iEvent.put(vector_mus_pid_TM2DCompatibilityTight	, "muspidTM2DCompatibilityTight" );
     iEvent.put(vector_mus_gfit_chi2			, "musgfitchi2"                  );
     iEvent.put(vector_mus_gfit_ndof			, "musgfitndof"                  );
     iEvent.put(vector_mus_gfit_validHits		, "musgfitvalidHits"             );
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonMaker);
