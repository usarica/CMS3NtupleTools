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
// $Id: PFCandidateMaker.cc,v 1.8 2011/09/14 17:41:57 slava77 Exp $
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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CMS2/NtupleMaker/interface/PFCandidateMaker.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

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

     pfCandidatesTag_		= iConfig.getParameter<InputTag>	("pfCandidatesTag");
     pfElectronsTag_		= iConfig.getParameter<InputTag>	("pfElectronsTag");
     tracksInputTag_            = iConfig.getParameter<InputTag>        ("tracksInputTag");
     minDR_electron_            = iConfig.getParameter<double>          ("minDRelectron");

     produces<vector<LorentzVector>	> ("pfcandsp4"              ).setBranchAlias("pfcands_p4"			          );
     produces<vector<LorentzVector>	> ("pfcandsposAtEcalp4"		  ).setBranchAlias("pfcands_posAtEcal_p4"		  );
     produces<vector<bool> >          ("pfcandsisMuIso"         ).setBranchAlias("pfcands_isMuIso"	        );
     produces<vector<float>	>         ("pfcandsecalE"	  	      ).setBranchAlias("pfcands_ecalE"		        );
     produces<vector<float>	>         ("pfcandshcalE"		        ).setBranchAlias("pfcands_hcalE"		        );
     produces<vector<float>	>         ("pfcandsrawEcalE"		    ).setBranchAlias("pfcands_rawEcalE"		      );
     produces<vector<float>	>         ("pfcandsrawHcalE"		    ).setBranchAlias("pfcands_rawHcalE"		      );
     produces<vector<float>	>         ("pfcandspS1E"		        ).setBranchAlias("pfcands_pS1E"			        );
     produces<vector<float>	>         ("pfcandspS2E"		        ).setBranchAlias("pfcands_pS2E"			        );
     produces<vector<float>	>         ("pfcandsdeltaP"		      ).setBranchAlias("pfcands_deltaP"		        );
     produces<vector<float>	>         ("pfcandsmvaepi"		      ).setBranchAlias("pfcands_mva_epi"		      );
     produces<vector<float>	>         ("pfcandsmvaemu"		      ).setBranchAlias("pfcands_mva_emu"		      );
     produces<vector<float>	>         ("pfcandsmvapimu"		      ).setBranchAlias("pfcands_mva_pimu"		      );
     produces<vector<float>	>         ("pfcandsmvanothinggamma"	).setBranchAlias("pfcands_mva_nothing_gamma");
     produces<vector<float>	>         ("pfcandsmvanothingnh"	  ).setBranchAlias("pfcands_mva_nothing_nh"	  );
     produces<vector<int>	>           ("pfcandscharge"		      ).setBranchAlias("pfcands_charge"		        );
     produces<vector<int> >           ("pfcandsparticleId"		  ).setBranchAlias("pfcands_particleId"		    );
     produces<vector<int>	>           ("pfcandsflag"		        ).setBranchAlias("pfcands_flag"			        );
     produces<vector<int>	>           ("pfcandstrkidx"		      ).setBranchAlias("pfcands_trkidx"		        );
     produces<vector<int>	>           ("pfcandspfmusidx"		    ).setBranchAlias("pfcands_pfmusidx"		      );
     produces<vector<int>	>           ("pfcandspfelsidx"		    ).setBranchAlias("pfcands_pfelsidx"		      );

     produces<float>                  ("evtfixgridrhoctr"       ).setBranchAlias("evt_fixgrid_rho_ctr"      );
     produces<float>                  ("evtfixgridrhofwd"       ).setBranchAlias("evt_fixgrid_rho_fwd"      );
     produces<float>                  ("evtfixgridrhoall"       ).setBranchAlias("evt_fixgrid_rho_all"      );
}

PFCandidateMaker::~PFCandidateMaker() {}
void  PFCandidateMaker::beginRun(edm::Run&, const edm::EventSetup& es) {}
void PFCandidateMaker::beginJob() {}
void PFCandidateMaker::endJob()   {}

// ------------ method called to produce the data  ------------
void PFCandidateMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

     auto_ptr<vector<LorentzVector> >	pfcands_p4		            (new vector<LorentzVector>  );
     auto_ptr<vector<LorentzVector> >	pfcands_posAtEcal_p4	    (new vector<LorentzVector>	);
     auto_ptr<vector<bool> >		      pfcands_isMuIso	          (new vector<bool> 		      );
     auto_ptr<vector<float> >		      pfcands_ecalE		          (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_hcalE		          (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_rawEcalE	        (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_rawHcalE	        (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_pS1E		          (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_pS2E		          (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_deltaP		        (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_mva_epi		        (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_mva_emu		        (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_mva_pimu	        (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_mva_nothing_gamma (new vector<float>		      );
     auto_ptr<vector<float> >		      pfcands_mva_nothing_nh	  (new vector<float>		      );
     auto_ptr<vector<int> >		        pfcands_charge		        (new vector<int>		        );
     auto_ptr<vector<int> >		        pfcands_particleId	      (new vector<int>		        );
     auto_ptr<vector<int> >	          pfcands_flag		          (new vector<int>        	  );
     auto_ptr<vector<int> >	          pfcands_trkidx		        (new vector<int>        	  );
     auto_ptr<vector<int> >	          pfcands_pfmusidx	        (new vector<int>        	  );    
     auto_ptr<vector<int> >	          pfcands_pfelsidx	        (new vector<int>        	  );
     auto_ptr<float >	                evt_fixgrid_rho_ctr 	    (new float           	      );
     auto_ptr<float >	                evt_fixgrid_rho_fwd 	    (new float           	      );
     auto_ptr<float >	                evt_fixgrid_rho_all 	    (new float           	      );
    
     //get pfcandidates
     Handle<PFCandidateCollection> pfCandidatesHandle;
     iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
     pfCandidates  = pfCandidatesHandle.product();

     //get pfelectrons
     typedef edm::ValueMap<reco::PFCandidatePtr> PFCandMap;
     Handle<PFCandMap> pfElectronsHandle;
     iEvent.getByLabel(pfElectronsTag_, pfElectronsHandle);
     const PFCandMap *pfElectrons  = pfElectronsHandle.product();
  

     // get tracks
     Handle<reco::TrackCollection>  track_h;
     iEvent.getByLabel(tracksInputTag_, track_h);

     int npfmus = 0;
  

     
       
      int iCand = 0;
      for( PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++ ) {
  
        //
	      int pfflags = 0;
	      for( unsigned int i = 0; i < 17; i++ ) {
          if(pf_it->flag((PFCandidate::Flags)i)) pfflags |= (1<<i);
	      }

        //
        bool isMuIso = false;
        const reco::PFCandidate &cand = (*pfCandidatesHandle)[iCand];
        if( cand.particleId() == 3 ) isMuIso = PFMuonAlgo::isIsolatedMuon( cand.muonRef() );
        iCand++;

  	    pfcands_p4			          ->push_back( LorentzVector(pf_it->px(), pf_it->py(), pf_it->pz(), pf_it->p())	  );
  	    pfcands_posAtEcal_p4		  ->push_back( LorentzVector(pf_it->positionAtECALEntrance().x(), pf_it->positionAtECALEntrance().y(), pf_it->positionAtECALEntrance().z(), 0.0) );
  	    pfcands_isMuIso		        ->push_back( isMuIso                                                            );
  	    pfcands_ecalE			        ->push_back( isfinite( pf_it->ecalEnergy() ) ? pf_it->ecalEnergy() : -9999.     );
  	    pfcands_hcalE			        ->push_back( pf_it->hcalEnergy()				                                        );
  	    pfcands_rawEcalE		      ->push_back( pf_it->rawEcalEnergy()			                                        );
  	    pfcands_rawHcalE		      ->push_back( pf_it->rawHcalEnergy()			                                        );
  	    pfcands_pS1E			        ->push_back( pf_it->pS1Energy()					                                        );
  	    pfcands_pS2E			        ->push_back( pf_it->pS2Energy()					                                        );
  	    pfcands_deltaP		        ->push_back( pf_it->deltaP()						                                        );
  	    pfcands_mva_epi		        ->push_back( pf_it->mva_e_pi()					                                        );
  	    pfcands_mva_emu		        ->push_back( pf_it->mva_e_mu()					                                        );
  	    pfcands_mva_pimu		      ->push_back( pf_it->mva_pi_mu()				                                          );
  	    pfcands_mva_nothing_gamma	->push_back( pf_it->mva_nothing_gamma()                                         );
  	    pfcands_mva_nothing_nh	  ->push_back( pf_it->mva_nothing_nh()		                                        );
  	    pfcands_charge		        ->push_back( pf_it->charge()						                                        );
  	    pfcands_particleId		    ->push_back( pf_it->translateTypeToPdgId(pf_it->particleId())	                  );
  	    pfcands_flag			        ->push_back( pfflags                                                            ); 


          //for charged pfcandidates, find corresponding track index
          //here we take the track directly from PFCandidate::trackRef()
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

          //find corresponding PFMuon index
          if (pf_it->particleId() == PFCandidate::mu){
            pfcands_pfmusidx->push_back(npfmus);
            ++npfmus;
          }else{
            pfcands_pfmusidx->push_back(-2);
          }

          //find corresponding PFElectron index
          if (pf_it->particleId() == PFCandidate::e){
	    int index = -1; 
	    if (pf_it->gsfTrackRef().isNonnull()) {
	      int pfGsfTkId = pf_it->gsfTrackRef().key();
	      unsigned int elsIndex = 0;
	      PFCandMap::const_iterator el_pit = pfElectrons->begin();
	      unsigned int nE = el_pit.size();
	      for(unsigned int iE = 0; iE < nE; ++iE){
		const PFCandidatePtr& el_it = el_pit[iE];
		if (el_it.isNull()) continue;
		int elGsfTkId = -1;
		if (el_it->gsfTrackRef().isNonnull()) elGsfTkId = el_it->gsfTrackRef().key();
		if (elGsfTkId==pfGsfTkId) {
		  index = elsIndex;
		  break;
		}
		elsIndex++;
	      }            
	    }
	    pfcands_pfelsidx->push_back(index);
	  } else {
            pfcands_pfelsidx->push_back(-2);
          }
 
          
     }//loop over candidate collection


     //define the phi bins
     vector<float> phibins;
     for (int i=0;i<10;i++) phibins.push_back(-TMath::Pi()+(2*i+1)*TMath::TwoPi()/20.);
     //define the eta bins
     vector<float> etabins_ctr;
     for (int i=0;i<8;++i) etabins_ctr.push_back(-2.1+0.6*i);
     vector<float> etabins_fwd;
     for (int i=0;i<10;++i) {
       if (i<5) etabins_fwd.push_back(-5.1+0.6*i);
       else etabins_fwd.push_back(2.7+0.6*(i-5));
     }
     vector<float> etabins_all;
     for (int i=0;i<18;++i) etabins_all.push_back(-5.1+0.6*i);
     //compute it
     *evt_fixgrid_rho_ctr = getFixGridRho(etabins_ctr,phibins);
     *evt_fixgrid_rho_fwd = getFixGridRho(etabins_fwd,phibins);
     *evt_fixgrid_rho_all = getFixGridRho(etabins_all,phibins);

     //
     iEvent.put(pfcands_p4,			            "pfcandsp4"		          );
     iEvent.put(pfcands_posAtEcal_p4,		    "pfcandsposAtEcalp4"	  );
     iEvent.put(pfcands_isMuIso,	          "pfcandsisMuIso"	      );
     iEvent.put(pfcands_ecalE,			        "pfcandsecalE"		      );
     iEvent.put(pfcands_hcalE,			        "pfcandshcalE"		      );
     iEvent.put(pfcands_rawEcalE,		        "pfcandsrawEcalE"	      );
     iEvent.put(pfcands_rawHcalE,		        "pfcandsrawHcalE"	      );
     iEvent.put(pfcands_pS1E,			          "pfcandspS1E"		        );
     iEvent.put(pfcands_pS2E,			          "pfcandspS2E"		        );
     iEvent.put(pfcands_deltaP,			        "pfcandsdeltaP"		      );
     iEvent.put(pfcands_mva_epi,		        "pfcandsmvaepi"		      );
     iEvent.put(pfcands_mva_emu,		        "pfcandsmvaemu"		      );
     iEvent.put(pfcands_mva_pimu,	  	      "pfcandsmvapimu"	      );
     iEvent.put(pfcands_mva_nothing_gamma,  "pfcandsmvanothinggamma");
     iEvent.put(pfcands_mva_nothing_nh,		  "pfcandsmvanothingnh"	  );
     iEvent.put(pfcands_charge,			        "pfcandscharge"		      );
     iEvent.put(pfcands_particleId,		      "pfcandsparticleId"	    );
     iEvent.put(pfcands_flag,			          "pfcandsflag"		        );
     iEvent.put(pfcands_trkidx,			        "pfcandstrkidx"		      );
     iEvent.put(pfcands_pfmusidx,		        "pfcandspfmusidx"	      );
     iEvent.put(pfcands_pfelsidx,		        "pfcandspfelsidx"	      );
     iEvent.put(evt_fixgrid_rho_ctr,		    "evtfixgridrhoctr"	    );
     iEvent.put(evt_fixgrid_rho_fwd,		    "evtfixgridrhofwd"	    );
     iEvent.put(evt_fixgrid_rho_all,		    "evtfixgridrhoall"	    );
 
}

float PFCandidateMaker::getFixGridRho(std::vector<float>& etabins,std::vector<float>& phibins) {
     float etadist = etabins[1]-etabins[0];
     float phidist = phibins[1]-phibins[0];
     float etahalfdist = (etabins[1]-etabins[0])/2.;
     float phihalfdist = (phibins[1]-phibins[0])/2.;
     vector<float> sumPFNallSMDQ;
     sumPFNallSMDQ.reserve(etabins.size()*phibins.size());
     for (unsigned int ieta=0;ieta<etabins.size();++ieta) {
       for (unsigned int iphi=0;iphi<phibins.size();++iphi) {
	 float pfniso_ieta_iphi = 0;
	 for(PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++) {
	   if (fabs(etabins[ieta]-pf_it->eta())>etahalfdist) continue;
	   if (fabs(reco::deltaPhi(phibins[iphi],pf_it->phi()))>phihalfdist) continue;
	   pfniso_ieta_iphi+=pf_it->pt();
	 }
	 sumPFNallSMDQ.push_back(pfniso_ieta_iphi);
       }
     }
     float evt_smdq = 0;
     sort(sumPFNallSMDQ.begin(),sumPFNallSMDQ.end());
     if (sumPFNallSMDQ.size()%2) evt_smdq = sumPFNallSMDQ[(sumPFNallSMDQ.size()-1)/2];
     else evt_smdq = (sumPFNallSMDQ[sumPFNallSMDQ.size()/2]+sumPFNallSMDQ[(sumPFNallSMDQ.size()-2)/2])/2.;
     return evt_smdq/(etadist*phidist);
}


//define this as a plug-in
DEFINE_FWK_MODULE(PFCandidateMaker);

