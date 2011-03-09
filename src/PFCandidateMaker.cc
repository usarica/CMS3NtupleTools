//-*- C++ -*-
//
// Package:    PFCandidateMaker
// Class:      PFCandidateMaker
// 
/**\class PFCandidateMaker PFCandidateMaker.cc CMS2/PFCandidateMaker/src/PFCandidateMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFCandidateMaker.cc,v 1.1 2011/03/09 08:53:11 benhoob Exp $
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
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CMS2/NtupleMaker/interface/PFCandidateMaker.h"


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
PFCandidateMaker::PFCandidateMaker(const edm::ParameterSet& iConfig) {

     pfCandidatesTag_		= iConfig.getParameter<InputTag>	("pfCandidatesTag"	);
     tracksInputTag_            = iConfig.getParameter<InputTag>        ("tracksInputTag");

     produces<vector<LorentzVector>	> ("pfcandsp4"                  ).setBranchAlias("pfcands_p4"			);
     produces<vector<LorentzVector>	> ("pfcandsposAtEcalp4"		).setBranchAlias("pfcands_posAtEcal_p4"		);
     produces<vector<float>		> ("pfcandsecalE"		).setBranchAlias("pfcands_ecalE"		);
     produces<vector<float>		> ("pfcandshcalE"		).setBranchAlias("pfcands_hcalE"		);
     produces<vector<float>		> ("pfcandsrawEcalE"		).setBranchAlias("pfcands_rawEcalE"		);
     produces<vector<float>		> ("pfcandsrawHcalE"		).setBranchAlias("pfcands_rawHcalE"		);
     produces<vector<float>		> ("pfcandspS1E"		).setBranchAlias("pfcands_pS1E"			);
     produces<vector<float>		> ("pfcandspS2E"		).setBranchAlias("pfcands_pS2E"			);
     produces<vector<float>		> ("pfcandsdeltaP"		).setBranchAlias("pfcands_deltaP"		);
     produces<vector<float>		> ("pfcandsmvaepi"		).setBranchAlias("pfcands_mva_epi"		);
     produces<vector<float>		> ("pfcandsmvaemu"		).setBranchAlias("pfcands_mva_emu"		);
     produces<vector<float>		> ("pfcandsmvapimu"		).setBranchAlias("pfcands_mva_pimu"		);
     produces<vector<float>		> ("pfcandsmvanothinggamma"	).setBranchAlias("pfcands_mva_nothing_gamma"	);
     produces<vector<float>		> ("pfcandsmvanothingnh"	).setBranchAlias("pfcands_mva_nothing_nh"	);
     produces<vector<int>		> ("pfcandscharge"		).setBranchAlias("pfcands_charge"		);
     produces<vector<int>		> ("pfcandsparticleId"		).setBranchAlias("pfcands_particleId"		);
     produces<vector<int>  		> ("pfcandsflag"		).setBranchAlias("pfcands_flag"			);
     produces<vector<int>  		> ("pfcandstrkidx"		).setBranchAlias("pfcands_trkidx"		);
}

PFCandidateMaker::~PFCandidateMaker() {

}


void  PFCandidateMaker::beginRun(edm::Run&, const edm::EventSetup& es) {
     

}

void PFCandidateMaker::beginJob() {


}

void PFCandidateMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void PFCandidateMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

     auto_ptr<vector<LorentzVector> >	pfcands_p4		   (new vector<LorentzVector>	);
     auto_ptr<vector<LorentzVector> >	pfcands_posAtEcal_p4	   (new vector<LorentzVector>	);
     auto_ptr<vector<float> >		pfcands_ecalE		   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_hcalE		   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_rawEcalE	   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_rawHcalE	   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_pS1E		   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_pS2E		   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_deltaP		   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_mva_epi		   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_mva_emu		   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_mva_pimu	   (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_mva_nothing_gamma  (new vector<float>		);
     auto_ptr<vector<float> >		pfcands_mva_nothing_nh	   (new vector<float>		);
     auto_ptr<vector<int> >		pfcands_charge		   (new vector<int>		);
     auto_ptr<vector<int> >		pfcands_particleId	   (new vector<int>		);
     auto_ptr<vector<int> >	        pfcands_flag		   (new vector<int>        	);
     auto_ptr<vector<int> >	        pfcands_trkidx		   (new vector<int>        	);
     auto_ptr<vector<float> >	        pfcands_isoChargedHadrons  (new vector<float>		);
     auto_ptr<vector<float> >	        pfcands_isoNeutralHadrons  (new vector<float>		);
     auto_ptr<vector<float> >	        pfcands_isoPhotons         (new vector<float>		);  

     //get pfcandidates
     Handle<PFCandidateCollection> pfCandidatesHandle;
     iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
     const PFCandidateCollection *pfCandidates  = pfCandidatesHandle.product();

     // get tracks
     Handle<reco::TrackCollection>  track_h;
     iEvent.getByLabel(tracksInputTag_, track_h);

     for(PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++)
     {
	  int pfflags = 0;
	  for(unsigned int i = 0; i < 17; i++)
	  {
	       if(pf_it->flag((PFCandidate::Flags)i))
		    pfflags |= (1<<i);
	  }

	  pfcands_p4			->push_back(LorentzVector(pf_it->px(), pf_it->py(), pf_it->pz(), pf_it->p())	);
	  pfcands_posAtEcal_p4		->push_back(LorentzVector(pf_it->positionAtECALEntrance().x(),
								  pf_it->positionAtECALEntrance().y(),
								  pf_it->positionAtECALEntrance().z(),
								  0.0)							);
	  pfcands_ecalE			->push_back(pf_it->ecalEnergy()							);
	  pfcands_hcalE			->push_back(pf_it->hcalEnergy()							);
	  pfcands_rawEcalE		->push_back(pf_it->rawEcalEnergy()						);
	  pfcands_rawHcalE		->push_back(pf_it->rawHcalEnergy()						);
	  pfcands_pS1E			->push_back(pf_it->pS1Energy()							);
	  pfcands_pS2E			->push_back(pf_it->pS2Energy()							);
	  pfcands_deltaP		->push_back(pf_it->deltaP()							);
	  pfcands_mva_epi		->push_back(pf_it->mva_e_pi()							);
	  pfcands_mva_emu		->push_back(pf_it->mva_e_mu()							);
	  pfcands_mva_pimu		->push_back(pf_it->mva_pi_mu()							);
	  pfcands_mva_nothing_gamma	->push_back(pf_it->mva_nothing_gamma()						);
	  pfcands_mva_nothing_nh	->push_back(pf_it->mva_nothing_nh()						);
	  pfcands_charge		->push_back(pf_it->charge()							);
	  pfcands_particleId		->push_back(pf_it->translateTypeToPdgId(pf_it->particleId())			);
	  pfcands_flag			->push_back(pfflags                                                             ); 


          //for charged pfcandidates, find corresponding track index
          if( pf_it->charge() != 0 ){
            
            reco::TrackRef pftrack = pf_it->trackRef();

            int trkidx = 0;
            bool foundTrack = false;

            reco::TrackCollection::const_iterator tracks_end = track_h->end();

            for (reco::TrackCollection::const_iterator itrk = track_h->begin(); itrk != tracks_end; ++itrk) {
              
              reco::TrackRef trkref( track_h , itrk - track_h->begin() );

              if( pftrack.key() == trkref.key() ){
              
                //sanity check
                float dpt = pftrack->pt() - trkref->pt();
                if( fabs( dpt ) > 0.1 ){
                  cout << "Warning: pfcandidate track pt - matched track pt = " << dpt << ", possible mismatch" << endl;
                }

                //found corresponding track
                pfcands_trkidx->push_back( trkidx );
                foundTrack = true;
                break;
              }
              
              ++trkidx;
            }

            if( !foundTrack ){
              //no matched track found, set trkidx to -1
              pfcands_trkidx->push_back(-1);
            }
            
          }else{
            //neutral particle, set trkidx to -2
            pfcands_trkidx->push_back(-2);
          }

     }//loop over candidate collection

     iEvent.put(pfcands_p4,			"pfcandsp4"		    );
     iEvent.put(pfcands_posAtEcal_p4,		"pfcandsposAtEcalp4"	    );
     iEvent.put(pfcands_ecalE,			"pfcandsecalE"		    );
     iEvent.put(pfcands_hcalE,			"pfcandshcalE"		    );
     iEvent.put(pfcands_rawEcalE,		"pfcandsrawEcalE"	    );
     iEvent.put(pfcands_rawHcalE,		"pfcandsrawHcalE"	    );
     iEvent.put(pfcands_pS1E,			"pfcandspS1E"		    );
     iEvent.put(pfcands_pS2E,			"pfcandspS2E"		    );
     iEvent.put(pfcands_deltaP,			"pfcandsdeltaP"		    );
     iEvent.put(pfcands_mva_epi,		"pfcandsmvaepi"		    );
     iEvent.put(pfcands_mva_emu,		"pfcandsmvaemu"		    );
     iEvent.put(pfcands_mva_pimu,		"pfcandsmvapimu"	    );
     iEvent.put(pfcands_mva_nothing_gamma,	"pfcandsmvanothinggamma"    );
     iEvent.put(pfcands_mva_nothing_nh,		"pfcandsmvanothingnh"	    );
     iEvent.put(pfcands_charge,			"pfcandscharge"		    );
     iEvent.put(pfcands_particleId,		"pfcandsparticleId"	    );
     iEvent.put(pfcands_flag,			"pfcandsflag"		    );
     iEvent.put(pfcands_trkidx,			"pfcandstrkidx"		    );
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFCandidateMaker);

