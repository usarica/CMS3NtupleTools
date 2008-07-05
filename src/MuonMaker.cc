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
// $Id: MuonMaker.cc,v 1.7 2008/07/05 23:04:49 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoMuon/MuonIdentification/interface/IdGlobalFunctions.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

MuonMaker::MuonMaker(const edm::ParameterSet& iConfig)
{
     // mu track quantities
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
     produces<vector<int> >             ("mustrkrefkey"         ).setBranchAlias("mus_trkrefkey"        );      // index of track from track ref stored in muon collection
     // muon quantities
     produces<vector<int> >	("musnmatches"			).setBranchAlias("mus_nmatches"                     	);	// number of stations with matched segments                                                  
     produces<vector<float> >	("museem"			).setBranchAlias("mus_e_em"                         	);	// energy in crossed ECAL crystalls                                                          
     produces<vector<float> >	("musehad"			).setBranchAlias("mus_e_had"                        	);	// energy in crossed HCAL towers                                                             
     produces<vector<float> >	("museho"			).setBranchAlias("mus_e_ho"                         	);	// energy in crossed HO towers                                                               
     produces<vector<float> >	("museemS9"			).setBranchAlias("mus_e_emS9"                       	);	// energy in 3x3 ECAL crystall shape                                                         
     produces<vector<float> >	("musehadS9"			).setBranchAlias("mus_e_hadS9"                      	);	// energy in 3x3 HCAL towers                                                                 
     produces<vector<float> >	("musehoS9"			).setBranchAlias("mus_e_hoS9"                       	);	// energy in 3x3 HO towers                                                                   
     produces<vector<float> >	("musiso03sumPt"		).setBranchAlias("mus_iso03_sumPt"                  	);	// sum of track Pt for cone of 0.3                                                           
     produces<vector<float> >	("musiso03emEt"			).setBranchAlias("mus_iso03_emEt"                   	);	// sum of ecal Et for cone of 0.3                                                            
     produces<vector<float> >	("musiso03hadEt"		).setBranchAlias("mus_iso03_hadEt"                  	);	// sum of hcal Et for cone of 0.3                                                            
     produces<vector<float> >	("musiso03hoEt"			).setBranchAlias("mus_iso03_hoEt"                   	);	// sum of ho Et for cone of 0.3                                                              
     produces<vector<int> >	("musiso03ntrk"			).setBranchAlias("mus_iso03_ntrk"                   	);	// number of tracks in the cone of 0.3                                                       
     produces<vector<float> >	("musiso05sumPt"		).setBranchAlias("mus_iso05_sumPt"                  	);	// sum of track Pt for cone of 0.5                                                           
     produces<vector<float> >	("musiso05emEt"			).setBranchAlias("mus_iso05_emEt"                   	);	// sum of ecal Et for cone of 0.5                                                            
     produces<vector<float> >	("musiso05hadEt"		).setBranchAlias("mus_iso05_hadEt"                  	);	// sum of hcal Et for cone of 0.5                                                            
     produces<vector<float> >	("musiso05hoEt"			).setBranchAlias("mus_iso05_hoEt"                   	);	// sum of ho Et for cone of 0.5                                                              
     produces<vector<int> >	("musiso05ntrk"			).setBranchAlias("mus_iso05_ntrk"                   	);	// number of tracks in the cone of 0.5                                                       
     produces<vector<int> >	("muspidTMLastStationLoose"	).setBranchAlias("mus_pid_TMLastStationLoose"       	);	// loose tracker muon identification based on muon/hadron penetration depth difference       
     produces<vector<int> >	("muspidTMLastStationTight"	).setBranchAlias("mus_pid_TMLastStationTight"       	);	// tight tracker muon identification based on muon/hadron penetration depth difference       
     produces<vector<int> >	("muspidTM2DCompatibilityLoose"	).setBranchAlias("mus_pid_TM2DCompatibilityLoose"   	);	// loose tracker muon likelihood identification based on muon matches and calo depositions   
     produces<vector<int> >	("muspidTM2DCompatibilityTight"	).setBranchAlias("mus_pid_TM2DCompatibilityTight"   	);	// tight tracker muon likelihood identification based on muon matches and calo depositions   
     produces<vector<float> >	("musgfitchi2"			).setBranchAlias("mus_gfit_chi2"                    	);	// chi2 of the global muon fit                                                               
     produces<vector<float> >	("musgfitndof"			).setBranchAlias("mus_gfit_ndof"                    	);	// number of degree of freedom of the global muon fit                                        
     produces<vector<int> >	("musgfitvalidHits"		).setBranchAlias("mus_gfit_validHits"               	);	// number of valid hits of the global muon fit                

     genParticlesInputTag = iConfig.getParameter<edm::InputTag>("genParticlesInputTag");                               
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
     std::auto_ptr<vector<int> >                vector_mus_trkrefkey    (new vector<int>                );
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

     // get MC particle collection
     edm::Handle<reco::GenParticleCollection> genParticlesHandle;
     iEvent.getByLabel(genParticlesInputTag, genParticlesHandle);
     
     for (edm::View<reco::Muon>::const_iterator muon = muon_h->begin(); 
	  muon != muons_end; ++muon) {

       //MC matching stuff
       const reco::GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*muon,
										    genParticlesHandle.product());
       int mcid = -999, mom_mcid = -999;
       LorentzVector mc_p4(0,0,0,0);
       if(matchedGenParticle != 0) {
	 mcid = matchedGenParticle->pdgId();
	 mc_p4 = matchedGenParticle->p4();
	 mom_mcid = MCUtilities::motherID(*matchedGenParticle)->pdgId();
       }

       // fill vectors
       vector_mus_p4           ->push_back(muon->p4());
       vector_mus_trk_p4       ->push_back(LorentzVector( muon->track().get()->px(), muon->track().get()->py(), muon->track().get()->pz(), muon->track().get()->p() ));
       vector_mus_d0           ->push_back(muon->track().isNonnull() ? muon->track()->d0()                       :	-999);
       vector_mus_z0           ->push_back(muon->track().isNonnull() ? muon->track()->dz()                       :	-999);
       vector_mus_vertexphi    ->push_back(muon->track().isNonnull() ? atan2( muon->track()->vy(), muon->track()->vx() )	:	-999);
       vector_mus_chi2         ->push_back(muon->track().isNonnull() ? muon->track()->chi2()                     :	-999);
       vector_mus_ndof         ->push_back(muon->track().isNonnull() ? muon->track()->ndof()                     :	-999);
       vector_mus_validHits    ->push_back(muon->track().isNonnull() ? muon->track()->numberOfValidHits()        :	-999);
       vector_mus_lostHits     ->push_back(muon->track().isNonnull() ? muon->track()->numberOfLostHits()         :	-999);
       vector_mus_d0Err        ->push_back(muon->track().isNonnull() ? muon->track()->d0Error()                  :	-999);
       vector_mus_z0Err        ->push_back(muon->track().isNonnull() ? muon->track()->dzError()                  :	-999);
       vector_mus_ptErr        ->push_back(muon->track().isNonnull() ? muon->track()->ptError()                  :	-999);
       vector_mus_etaErr       ->push_back(muon->track().isNonnull() ? muon->track()->etaError()                 :	-999);
       vector_mus_phiErr       ->push_back(muon->track().isNonnull() ? muon->track()->phiError()                 :	-999);
       vector_mus_mc_p4        ->push_back(mc_p4                                                                            );
       vector_mus_mc_id        ->push_back(mcid                                                                             );
       vector_mus_charge       ->push_back(muon->charge()                                                                   );
       vector_mus_mc_motherid  ->push_back(mom_mcid                                                                         );
       vector_mus_outerPhi     ->push_back(-999);
       vector_mus_outerEta     ->push_back(-999);
       vector_mus_trkrefkey    ->push_back(muon->track().index());
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
     iEvent.put(vector_mus_trkrefkey    , "mustrkrefkey"         );
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
